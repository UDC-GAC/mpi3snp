/*
 * GPUEngine.cu
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#include "GPUEngine.h"
#include "ThreadParams.h"
#include "EntropySearch.h"
#include <float.h>

GPUEngine::GPUEngine(Options *options) {
    /*check the availability of GPUs*/
    if (GPUInfo::getGPUInfo()->getNumGPUs() == 0) {
        Utils::exit("No compatible GPUs are available in your machine\n");
    }

//    if (options->isHeteroGPUs()) {
//        distributor = new GPUSNPDistributorStatic(options);
//    } else {
//        distributor = new GPUSNPDistributor(options);
//    }
    distributor = new GPUSNPDistributor(options);
    gpu_ids = options->Get_GPU_Ids();
    is_mi = options->isMI();
    num_outputs = options->getNumOutputs();
}

GPUEngine::~GPUEngine() {
    delete distributor;
}

void GPUEngine::run(std::vector<MutualInfo> &mutual_info, Statistics &statistics) {
    std::string snp_load_label("SNPs load time");
    statistics.Begin_timer(snp_load_label);
    distributor->loadSNPSet();
    statistics.End_timer(snp_load_label);

    vector<pthread_t> threadIDs(gpu_ids.size(), 0);
    vector<ThreadParams *> threadParams(gpu_ids.size());
    for (int tid = 0; tid < gpu_ids.size(); tid++) {
        // Create parameters for CPU threads
        threadParams[tid] = new ThreadParams(tid, num_outputs, distributor, gpu_ids[tid], is_mi, statistics);
        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, handle, threadParams[tid]) != 0) {
            Utils::exit("Thread creating failed\n");
        }
    }

    MutualInfo *auxMutualInfo = new MutualInfo[gpu_ids.size() * num_outputs];

    // Wait for the completion of all threads
    for (int tid = 0; tid < gpu_ids.size(); tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * num_outputs], threadParams[tid]->_mutualInfo,
               num_outputs * sizeof(MutualInfo));
        delete threadParams[tid];
    }

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo, auxMutualInfo + num_outputs * gpu_ids.size());
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], auxMutualInfo + num_outputs * (gpu_ids.size() - 1), sizeof(MutualInfo) * num_outputs);


#ifdef DEBUG
    uint32_t numAnalyzed = 0;

    for(int tid=0; tid<_options->getNumGPUs(); tid++){
        numAnalyzed+=_threadParams[tid]->_numAnalyzed;
    }

    Utils::log("Total analysis: %" PRIu64 "\n", numAnalyzed);
#endif

    // Release the distributor
    delete[] auxMutualInfo;
}

void *GPUEngine::handle(void *arg) {
    ThreadParams *params = (ThreadParams *) arg;
    GPUSNPDistributor *distributor = params->_distributor;
    uint16_t num_outputs = params->_numOutputs;
    int gpu_id = params->_gpu;
    bool isMI = params->_isMI;
    Statistics &statistics = params->statistics;

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
    std::string timer_label;
    timer_label += "GPU " + std::to_string(gpu_id) + " runtime";
    std::string analysis_label;
    analysis_label += "GPU " + std::to_string(gpu_id) + " analysis";

    statistics.Begin_timer(timer_label);
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

    statistics.End_timer(timer_label);
    statistics.Add_value(analysis_label, myTotalAnal);
#endif

    cudaFreeHost(auxIds);
    myCheckCudaError;
    delete search;
    return NULL;
}