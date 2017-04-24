/*
 * SearchMI.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#include "SearchMI.h"

SearchMI::SearchMI(Options *options) : Search(options) {
    _engines.resize(options->getNumThreads());
    _threadParams.resize(options->getNumThreads());

    for (int tid = 0; tid < options->getNumThreads(); ++tid) {
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
    double stime = Utils::getSysTime();
    double etime;

    int proc_id, num_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    // Declaration of MPI Datatypes
    MPI_Datatype MPI_MUTUAL_INFO;
    MPI_Type_contiguous(sizeof(MutualInfo), MPI_CHAR, &MPI_MUTUAL_INFO);
    MPI_Type_commit(&MPI_MUTUAL_INFO);

    SNPDistributor *distributor = new SNPDistributor(_options);
    distributor->loadSNPSet();

    etime = Utils::getSysTime();
    IOMpi::Instance().Cprintf("Loaded %ld SNPs in %.2f seconds\n", distributor->getNumSnp(), etime - stime);

    vector<pthread_t> threadIDs(_options->getNumThreads(), 0);

    // Computation of the single-SNP entropy
    for (int tid = 0; tid < _options->getNumThreads(); tid++) {
        // All threads use the same distributor so it gives them the correct pair of SNPs
        _threadParams[tid]->init(distributor, _options->getNumOutputs());

        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, _threadMI, _threadParams[tid]) != 0) {
            Utils::exit("Thread creating failed\n");
        }
    }

    MutualInfo *auxMutualInfo;
    if (proc_id == 0) { // Master
        auxMutualInfo = new MutualInfo[num_proc * _options->getNumThreads() * _options->getNumOutputs()];
    } else {
        auxMutualInfo = new MutualInfo[_options->getNumThreads() * _options->getNumOutputs()];
    }

    // Wait for the completion of all threads
    for (int tid = 0; tid < _options->getNumThreads(); tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * _options->getNumOutputs()], _threadParams[tid]->_mutualInfo,
               _options->getNumOutputs() * sizeof(MutualInfo));
    }

    // Gather all the results on the master process and print output
    if (proc_id == 0) {
        for (int rank = 1; rank < num_proc; rank++) {
            MPI_Recv(&auxMutualInfo[rank * _options->getNumThreads() * _options->getNumOutputs()],
                     _options->getNumThreads() * _options->getNumOutputs(), MPI_MUTUAL_INFO, rank, 123,
                     MPI_COMM_WORLD,
                     NULL);
        }

        // Sort the auxiliar array and print the results
        std::sort(auxMutualInfo, auxMutualInfo + _options->getNumOutputs() * _options->getNumThreads() * num_proc);
        distributor->printMI(auxMutualInfo + _options->getNumOutputs() * (_options->getNumThreads() * num_proc - 1),
                             _options->getNumOutputs());
    } else {
        MPI_Send(auxMutualInfo, _options->getNumThreads() * _options->getNumOutputs(), MPI_MUTUAL_INFO, 0, 123,
                 MPI_COMM_WORLD);
    }

#ifdef BENCHMARKING
    uint32_t numAnalyzed = 0;

    for(int tid=0; tid<_options->getNumThreads(); tid++){
        IOMpi::Instance().Cprintf("CPU thread (%d) %f seconds calculating %lu analysis\n",
                tid, _threadParams[tid]->_runtime, _threadParams[tid]->_numAnalyzed);
        numAnalyzed+=_threadParams[tid]->_numAnalyzed;
    }

    IOMpi::Instance().Cprintf("Total analysis: %" PRIu64 "\n", numAnalyzed);
#endif

    // Release the distributor
    MPI_Type_free(&MPI_MUTUAL_INFO);
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
    double stime = Utils::getSysTime();

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
    params->_runtime = Utils::getSysTime() - stime;
    memcpy(params->_mutualInfo, mutualInfo, numOutputs * sizeof(MutualInfo));

    delete[] auxIds;
    delete[] mutualInfo;
    return NULL;
}
