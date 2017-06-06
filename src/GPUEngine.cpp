/*
 * GPUEngine.cu
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#include "GPUEngine.h"
#include "GPUSNPDistributorStatic.h"
#include "ThreadParams.h"
#include "EntropySearch.h"
#include <float.h>

GPUEngine::GPUEngine(Options *options) {
    if (options->isHeteroGPUs()) {
        distributor = new GPUSNPDistributorStatic(options);
    } else {
        distributor = new GPUSNPDistributor(options);
    }
    num_gpus = options->getNumGPUs();
    gpu_id = new int[num_gpus];
    for (int i=0; i<num_gpus; i++){
        gpu_id[i] = options->getGPUId(i);
    }
    is_mi = options->isMI();
    num_outputs = options->getNumOutputs();
}

GPUEngine::~GPUEngine() {
    delete distributor;
    delete[] gpu_id;
}

void GPUEngine::run() {
    double stime = Utils::getSysTime();
    double etime;

    distributor->loadSNPSet();

    etime = Utils::getSysTime();
    Utils::log("Loaded %ld SNPs in %.2f seconds\n", distributor->getNumSnp(), etime - stime);

    vector<pthread_t> threadIDs(num_gpus, 0);
    vector<ThreadParams *> threadParams(num_gpus);
    for (int tid = 0; tid < num_gpus; tid++) {
        // Create parameters for CPU threads
        threadParams[tid] = new ThreadParams(tid, num_outputs, distributor, gpu_id[tid], is_mi);
        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, handle, threadParams[tid]) != 0) {
            Utils::exit("Thread creating failed\n");
        }
    }

    MutualInfo *auxMutualInfo = new MutualInfo[num_gpus * num_outputs];

    // Wait for the completion of all threads
    for (int tid = 0; tid < num_gpus; tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * num_outputs], threadParams[tid]->_mutualInfo,
               num_outputs * sizeof(MutualInfo));
        delete threadParams[tid];
    }

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo, auxMutualInfo + num_outputs * num_gpus);
    distributor->printMI(auxMutualInfo + num_outputs * (num_gpus - 1), num_outputs);

    Utils::log("3-SNP analysis finalized\n");

#ifdef DEBUG
    uint32_t numAnalyzed = 0;

    for(int tid=0; tid<_options->getNumGPUs(); tid++){
        numAnalyzed+=_threadParams[tid]->_numAnalyzed;
    }

    Utils::log("Total analysis: %" PRIu64 "\n", numAnalyzed);
#endif

    // Release the distributor
    delete distributor;
    delete[] auxMutualInfo;
}

void *GPUEngine::handle(void *arg) {
    ThreadParams *params = (ThreadParams *) arg;
    GPUSNPDistributor *distributor = params->_distributor;
    uint16_t num_outputs = params->_numOutputs;
    int gpu_id = params->_gpu;
    bool isMI = params->_isMI;

    GPUInfo::getGPUInfo()->setDevice(gpu_id);

    EntropySearch *search = new EntropySearch(isMI, distributor->getNumSnp(), distributor->getNumCases(),
                                                distributor->getNumCtrls(), num_outputs,
                                                distributor->getHost0Cases(), distributor->getHost1Cases(),
                                                distributor->getHost2Cases(), distributor->getHost0Ctrls(),
                                                distributor->getHost1Ctrls(), distributor->getHost2Ctrls());

    uint2 *auxIds;
    cudaMallocHost(&auxIds, NUM_PAIRS_BLOCK * sizeof(uint2));
    myCheckCudaError;

    // Variables to work with the outputs
    MutualInfo *mutualInfo = params->_mutualInfo;
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
        numPairsBlock = distributor->getPairsSNPs(auxIds, myTotalAnal, params->_tid);

        if (numPairsBlock) {
            search->mutualInfo(numPairsBlock, auxIds, mutualInfo, minMI, minMIPos, numEntriesWithMI);
        } else {
            moreAnal = false;
        }
    }

#ifdef BENCHMARKING
    cudaDeviceSynchronize();
    etime = Utils::getSysTime();
    Utils::log("GPU thread (%d) %f seconds calculating %lu analysis\n",
               params->_tid, etime - stime, myTotalAnal);
#endif

    params->_numAnalyzed = myTotalAnal;

    cudaFreeHost(auxIds);
    myCheckCudaError;
    delete search;
    return NULL;
}