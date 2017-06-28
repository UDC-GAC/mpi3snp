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
#include <cstring>

GPUEngine::GPUEngine(unsigned int proc_num, unsigned int proc_id, std::vector<unsigned int> gpu_ids, bool use_mi) :
        proc_num(proc_num),
        proc_id(proc_id),
        gpu_ids(gpu_ids),
        use_mi(use_mi) {
    int avail_gpus = 0;
    if (cudaSuccess != cudaGetDeviceCount(&avail_gpus))
        throw CUDAError(cudaGetLastError());
    if (avail_gpus == 0) {
        throw CUDAError("Could not find any CUDA-enabled GPU");
    }

    for (unsigned int id : gpu_ids) {
        cudaDeviceProp gpu_prop;
        if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, id))
            throw CUDAError(cudaGetLastError());
        if (gpu_prop.major < 2 || !gpu_prop.canMapHostMemory) {
            throw CUDAError("GPU " + std::to_string(id) + " does not meet compute capabilities\n" +
                            "Name: " + gpu_prop.name + "\n" + "Compute capability: " +
                            std::to_string(gpu_prop.major) + "." + std::to_string(gpu_prop.minor));
        }
    }
}

void GPUEngine::run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, unsigned int num_outputs,
                    Statistics &statistics) {
    std::string snp_load_label("SNPs load time");
    statistics.Begin_timer(snp_load_label);
    Dataset dataset(tped, tfam, Dataset::Transposed);
    statistics.End_timer(snp_load_label);

    Distributor distributor(proc_num, proc_id, dataset.Get_SNP_count(), gpu_ids.size() > 1);
    std::vector<pthread_t> threadIDs(gpu_ids.size(), 0);
    std::vector<ThreadParams *> threadParams(gpu_ids.size());
    for (int tid = 0; tid < gpu_ids.size(); tid++) {
        // Create parameters for CPU threads
        threadParams[tid] = new ThreadParams(gpu_ids[tid], num_outputs, dataset, distributor, use_mi, statistics);
        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, handle, threadParams[tid]) != 0) {
            //Utils::exit("Thread creating failed\n");
            // TODO: thread error handling
            exit(0);
        }
    }

    MutualInfo *auxMutualInfo = new MutualInfo[gpu_ids.size() * num_outputs];

    // Wait for the completion of all threads
    for (int tid = 0; tid < gpu_ids.size(); tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * num_outputs], threadParams[tid]->mutual_info,
               num_outputs * sizeof(MutualInfo));
        delete threadParams[tid];
    }

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo, auxMutualInfo + num_outputs * gpu_ids.size());
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], auxMutualInfo + num_outputs * (gpu_ids.size() - 1), sizeof(MutualInfo) * num_outputs);

    // Release the distributor
    delete[] auxMutualInfo;
}

void *GPUEngine::handle(void *arg) {
    ThreadParams *params = (ThreadParams *) arg;
    Dataset &dataset = params->dataset;
    Distributor &distributor = params->distributor;
    unsigned int num_outputs = params->num_outputs;
    unsigned int gpu_id = params->gpu_id;
    bool isMI = params->mi;
    Statistics &statistics = params->statistics;

    if (cudaSuccess != cudaSetDevice(gpu_id))
        throw CUDAError(cudaGetLastError());

    uint32_t *hCa0, *hCa1, *hCa2, *hCt0, *hCt1, *hCt2;
    if (cudaSuccess != cudaMallocHost(&hCa0, dataset.Get_cases()[0][0].size() * sizeof(uint32_t)))
        throw CUDAError(cudaGetLastError());
    memcpy(hCa0, &dataset.Get_cases()[0][0][0], dataset.Get_cases()[0][0].size() * sizeof(uint32_t));
    if (cudaSuccess != cudaMallocHost(&hCa1, dataset.Get_cases()[0][1].size() * sizeof(uint32_t)))
        throw CUDAError(cudaGetLastError());
    memcpy(hCa1, &dataset.Get_cases()[0][1][0], dataset.Get_cases()[0][1].size() * sizeof(uint32_t));
    if (cudaSuccess != cudaMallocHost(&hCa2, dataset.Get_cases()[0][2].size() * sizeof(uint32_t)))
        throw CUDAError(cudaGetLastError());
    memcpy(hCa2, &dataset.Get_cases()[0][2][0], dataset.Get_cases()[0][2].size() * sizeof(uint32_t));
    if (cudaSuccess != cudaMallocHost(&hCt0, dataset.Get_ctrls()[0][0].size() * sizeof(uint32_t)))
        throw CUDAError(cudaGetLastError());
    memcpy(hCt0, &dataset.Get_ctrls()[0][0][0], dataset.Get_ctrls()[0][0].size() * sizeof(uint32_t));
    if (cudaSuccess != cudaMallocHost(&hCt1, dataset.Get_ctrls()[0][1].size() * sizeof(uint32_t)))
        throw CUDAError(cudaGetLastError());
    memcpy(hCt1, &dataset.Get_ctrls()[0][1][0], dataset.Get_ctrls()[0][1].size() * sizeof(uint32_t));
    if (cudaSuccess != cudaMallocHost(&hCt2, dataset.Get_ctrls()[0][2].size() * sizeof(uint32_t)))
        throw CUDAError(cudaGetLastError());
    memcpy(hCt2, &dataset.Get_ctrls()[0][2][0], dataset.Get_ctrls()[0][2].size() * sizeof(uint32_t));

    EntropySearch *search = new EntropySearch(isMI, dataset.Get_SNP_count(), dataset.Get_case_count(),
                                              dataset.Get_ctrl_count(), num_outputs,
                                              hCa0, hCa1, hCa2, hCt0, hCt1, hCt2);

    uint2 *auxIds;
    if (cudaSuccess != cudaMallocHost(&auxIds, Distributor::DEFAULT_PAIRS_BLOCK * sizeof(uint2)))
        throw CUDAError(cudaGetLastError());

    // Variables to work with the outputs
    MutualInfo *mutualInfo = params->mutual_info;
    // The minimum value in the array
    float minMI = FLT_MAX;
    // The position of the minimum value
    uint16_t minMIPos = 0;
    // Number of entries of the array full
    uint16_t numEntriesWithMI = 0;

    bool moreAnal = true;
    uint64_t myTotalAnal = 0;
    uint64_t numPairsBlock = 0;

    std::string timer_label;
    timer_label += "GPU " + std::to_string(gpu_id) + " runtime";
    std::string analysis_label;
    analysis_label += "GPU " + std::to_string(gpu_id) + " analysis";

    statistics.Begin_timer(timer_label);

    while (moreAnal) {
        // Take some SNPs
        numPairsBlock = distributor.Get_pairs(auxIds, myTotalAnal);

        if (numPairsBlock) {
            search->mutualInfo(numPairsBlock, auxIds, mutualInfo, minMI, minMIPos, numEntriesWithMI);
        } else {
            moreAnal = false;
        }
    }

    cudaDeviceSynchronize();

    statistics.End_timer(timer_label);
    statistics.Add_value(analysis_label, myTotalAnal);

    if (cudaSuccess != cudaFreeHost(auxIds))
        throw CUDAError(cudaGetLastError());
    if (cudaSuccess != cudaFreeHost(hCa0))
        throw CUDAError(cudaGetLastError());
    if (cudaSuccess != cudaFreeHost(hCa1))
        throw CUDAError(cudaGetLastError());
    if (cudaSuccess != cudaFreeHost(hCa2))
        throw CUDAError(cudaGetLastError());
    if (cudaSuccess != cudaFreeHost(hCt0))
        throw CUDAError(cudaGetLastError());
    if (cudaSuccess != cudaFreeHost(hCt1))
        throw CUDAError(cudaGetLastError());
    if (cudaSuccess != cudaFreeHost(hCt2))
        throw CUDAError(cudaGetLastError());
    delete search;
    return nullptr;
}