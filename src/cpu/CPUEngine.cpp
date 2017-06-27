/*
 * Engine.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#include <cfloat>
#include "CPUEngine.h"
#include "ThreadParams.h"
#include "EntropySearch.h"

CPUEngine::CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi) {
    this->num_proc = num_proc;
    this->proc_id = proc_id;
    this->num_threads = num_threads;
    this->use_mi = use_mi;
}

CPUEngine::~CPUEngine() {

}

void CPUEngine::execute(std::string tped_file, std::string tfam_file, std::vector<MutualInfo> &mutual_info,
                        uint16_t num_outputs, Statistics &statistics) {
    double stime = Utils::getSysTime();
    double etime;

    SNPDistributor distributor(num_proc, proc_id, num_threads, tped_file, tfam_file);
    distributor.loadSNPSet();

    etime = Utils::getSysTime();
    IOMpi::Instance().Cprintf<IOMpi::D>("Loaded %ld SNPs (%ld/%ld cases/controls) in %.2f seconds\n",
                              distributor.getNumSnp(), distributor.getNumCases(), distributor.getNumCtrls(),
                              etime - stime);

    vector<pthread_t> threadIDs(num_threads, 0);

    // Computation of the single-SNP entropy
    std::vector<ThreadParams *> params(num_threads);
    for (int tid = 0; tid < num_threads; tid++) {
        params[tid] = new ThreadParams(tid);
        // All threads use the same distributor so it gives them the correct pair of SNPs
        params[tid]->init(&distributor, num_outputs);

        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, threadMI, params[tid]) != 0) {
            Utils::exit("Thread creating failed\n");
        }
    }

    MutualInfo *auxMutualInfo = new MutualInfo[num_threads * num_outputs];

    // Wait for the completion of all threads
    for (int tid = 0; tid < num_threads; tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * num_outputs], params[tid]->_mutualInfo,
               num_outputs * sizeof(MutualInfo));
        delete params[tid];
    }

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo, auxMutualInfo + num_outputs * num_threads);
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], auxMutualInfo + num_outputs * (num_threads - 1), sizeof(MutualInfo) * num_outputs);
}

void* CPUEngine::threadMI(void *arg) {
    ThreadParams *params = (ThreadParams *) arg;
    uint16_t numOutputs = params->_numOutputs;
    SNPDistributor *distributor = params->_distributor;

    EntropySearch search(distributor->getNumSnp(), distributor->getNumCases(), distributor->getNumCtrls());
    const vector<SNP2 *> &snpSet = distributor->getSnpSet();
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
            myTotalAnal += search.mutualInfo(snpSet, numPairsBlock, auxIds, mutualInfo, numOutputs, minMI,
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

