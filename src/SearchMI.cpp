/*
 * SearchMI.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 *
 *  Modified on: December 29, 2016
 *           by: Christian Ponte Fernández
 */

#include "SearchMI.h"

SearchMI::SearchMI(Options *options) : Search(options) {
    _engines.resize(options->getNumCPUs());
    _threadParams.resize(options->getNumCPUs());

    for (int tid = 0; tid < options->getNumCPUs(); ++tid) {
        _engines[tid] = new Engine(options);
        // Create parameters for CPU threads
        _threadParams[tid] = new ThreadParams(tid, _engines[tid]);
    }
}

SearchMI::~SearchMI() {
    for (int i = 0; i < _engines.size(); i++) {
        delete _engines[i];
        delete _threadParams[i];
    }

    _engines.clear();
    _threadParams.clear();
}

void SearchMI::execute() {
    // Obtain the whole SNPSet on all processes
    SNPDistributor *distributor = new SNPDistributor(_options);
    vector<SNP *> snpSet = distributor->getSnpSet();
    SNPDistributor::ClassSet_t classSet = distributor->getClassSet();
    delete distributor;

    MPI_Init(NULL, NULL);

    _mpiMI(_options, _threadParams, snpSet, classSet);

    MPI_Finalize();
}

void *SearchMI::_mpiMI(Options *options, vector<ThreadParams *> threadParams, vector<SNP *> snpSet,
                       SNPDistributor::ClassSet_t classSet) {
    double stime = Utils::getSysTime();
    double etime;

    int MPIRank, MPISize;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIRank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPISize);

    // Calculate lower and upper SNP index for each process
    uint32_t tuplesPerProc, n, k, sum;
    n = snpSet.size();
    tuplesPerProc = (n - 1) * n / (2 * MPISize);
    sum = 0;
    k = n;
    for (int procId = 0; procId <= MPIRank; procId++) {
        sum += n - k;
        n = k;
        k = (1 + sqrt(1 + 4 * (n * n - n - 2 * tuplesPerProc))) / 2;
    }
    if (MPIRank == MPISize - 1) {
        k = 0;
    }

    // Create one SNPDistributor per process
    SNPDistributor *distributor = new SNPDistributor(options, snpSet, classSet, sum, sum + n - k);

    etime = Utils::getSysTime();
    Utils::log("Process %i: loaded %ld SNPs in %.2f seconds, computing %i SNPs\n", MPIRank, snpSet.size(),
               etime - stime,
               n - k);

    vector<pthread_t> threadIDs(options->getNumCPUs(), 0);

    // Computation of the single-SNP entropy
    for (int tid = 0; tid < options->getNumCPUs(); tid++) {
        // All threads use the same distributor so it gives them the correct pair of SNPs
        threadParams[tid]->init(distributor, options->getNumOutputs());

        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, _threadMI, threadParams[tid]) != 0) {
            Utils::exit("Thread creating failed\n");
        }
    }

    MutualInfo *auxMutualInfo;
    if (MPIRank == 0) { // Master
        auxMutualInfo = new MutualInfo[MPISize * options->getNumCPUs() * options->getNumOutputs()];
    } else {
        auxMutualInfo = new MutualInfo[options->getNumCPUs() * options->getNumOutputs()];
    }

    // Wait for the completion of all threads
    for (int tid = 0; tid < options->getNumCPUs(); tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * options->getNumOutputs()], threadParams[tid]->_mutualInfo,
               options->getNumOutputs() * sizeof(MutualInfo));
    }

    // Gather all the results on the master process and print output
    if (MPIRank == 0) {
        for (int rank = 1; rank <= MPISize; rank++) {
            MPI_Recv(&auxMutualInfo[rank * options->getNumCPUs() * options->getNumOutputs()],
                     options->getNumCPUs() * options->getNumOutputs(), MPI_INT, rank, MPI_TAG_OUTPUT, MPI_COMM_WORLD,
                     NULL);
        }

        // Sort the auxiliar array and print the results
        std::sort(auxMutualInfo, auxMutualInfo + options->getNumOutputs() * options->getNumCPUs());
        distributor->printMI(auxMutualInfo + options->getNumOutputs() * (options->getNumCPUs() - 1),
                             options->getNumOutputs());
    } else {
        MPI_Send(auxMutualInfo, options->getNumCPUs() * options->getNumOutputs(), MPI_INT, 0, MPI_TAG_OUTPUT,
                 MPI_COMM_WORLD);
    }

    Utils::log("Process %i: 3-SNP analysis finalized\n", MPIRank);

#ifdef DEBUG
    uint32_t numAnalyzed = 0;

    for(int tid=0; tid<_options->getNumCPUs(); tid++){
        numAnalyzed+=_threadParams[tid]->_numAnalyzed;
    }

    Utils::log("Total analysis: %" PRIu64 "\n", numAnalyzed);
#endif

    // Release the distributor
    delete distributor;
    delete[] auxMutualInfo;
    threadIDs.clear();
}

void *SearchMI::_threadMI(void *arg) {
    ThreadParams *params = (ThreadParams *) arg;

    Engine *engine = params->_engine;

    SNPDistributor *distributor = params->_distributor;
    uint32_t numSnp = distributor->getNumSnp();
    uint16_t numCases = distributor->getNumCases();
    uint16_t numCtrls = distributor->getNumCtrls();
    uint16_t numOutputs = params->_numOutputs;
    engine->setNums(numSnp, numCases, numCtrls);

    vector<SNP *> snpSet = distributor->getSnpSet();
    // In this case the ids are necessary to access to the single-SNP entropy
    uint32_t *auxIds = new uint32_t[2 * NUM_PAIRS_BLOCK];

    // Variables to work with the outputs
    MutualInfo *mutualInfo = new MutualInfo[numOutputs];
    // The minimum value in the array
    float minMI = FLT_MAX;
    // The position of the minimum value
    uint16_t minMIPos = 0;
    // Number of entries of the array full
    uint16_t numEntriesWithMI = 0;

    bool moreAnal = true;
    uint64_t myTotalAnal = 0;

    uint64_t numPairsBlock = 0;

#ifdef BENCHMARKING
    double stime = Utils::getSysTime();
    double etime;
#endif

    while (moreAnal) {
        // Take some SNPs
        numPairsBlock = distributor->getPairsSNPs(auxIds);

        if (numPairsBlock) {
            myTotalAnal += engine->mutualInfo(snpSet, numPairsBlock, auxIds, mutualInfo, numOutputs, minMI,
                                              minMIPos, numEntriesWithMI);
        } else {
            moreAnal = false;
        }
    }

    params->_numAnalyzed = myTotalAnal;
    memcpy(params->_mutualInfo, mutualInfo, numOutputs * sizeof(MutualInfo));

#ifdef BENCHMARKING
    etime = Utils::getSysTime();
    Utils::log("CPU thread (%d) %f seconds calculating %lu analysis\n",
               params->_tid, etime - stime, myTotalAnal);
#endif

    delete[] auxIds;
    delete[] mutualInfo;
    return NULL;
}
