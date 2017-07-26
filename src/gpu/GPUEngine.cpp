/*
 * GPUEngine.cu
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#include "GPUEngine.h"
#include "EntropySearch.h"
#include "../Dataset.h"
#include "../Distributor.h"
#include <cfloat>
#include <cstring>

GPUEngine::GPUEngine(unsigned int proc_num, unsigned int proc_id, bool use_mi) :
        proc_num(proc_num),
        proc_id(proc_id),
        use_mi(use_mi) {
    int avail_gpus = 0;
    if (cudaSuccess != cudaGetDeviceCount(&avail_gpus))
        throw CUDAError(cudaGetLastError());
    if (avail_gpus == 0) {
        throw CUDAError("Could not find any CUDA-enabled GPU");
    }

    gpu_id = proc_id % avail_gpus;
    if (cudaSuccess != cudaSetDevice(gpu_id))
        throw CUDAError(cudaGetLastError());
    cudaDeviceProp gpu_prop;
    if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gpu_id))
        throw CUDAError(cudaGetLastError());
    if (gpu_prop.major < 2 || !gpu_prop.canMapHostMemory) {
        throw CUDAError("GPU " + std::to_string(gpu_id) + " does not meet compute capabilities\n" +
                        "Name: " + gpu_prop.name + "\n" + "Compute capability: " +
                        std::to_string(gpu_prop.major) + "." + std::to_string(gpu_prop.minor));
    }
}

void GPUEngine::run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, unsigned int num_outputs,
                    Statistics &statistics) {
    std::string snp_load_label("SNPs load time");
    statistics.Begin_timer(snp_load_label);
    Dataset dataset(tped, tfam, Dataset::Transposed);
    statistics.End_timer(snp_load_label);

    Distributor distributor(proc_num, proc_id, dataset.Get_SNP_count());

    const unsigned int pairs_iter = 5000;
    EntropySearch search(use_mi, dataset.Get_SNP_count(), dataset.Get_case_count(), dataset.Get_ctrl_count(),
                         num_outputs, pairs_iter, dataset.Get_cases(), dataset.Get_ctrls());
    uint2 *auxIds;
    if (cudaSuccess != cudaMallocHost(&auxIds, pairs_iter * sizeof(uint2)))
        throw CUDAError(cudaGetLastError());

    // The minimum value in the array
    float minMI = FLT_MAX;
    // The position of the minimum value
    uint16_t minMIPos = 0;
    // Number of entries of the array full
    uint16_t numEntriesWithMI = 0;

    uint64_t myTotalAnal = 0;
    uint64_t numPairsBlock = 0;

    std::string timer_label;
    timer_label += "GPU " + std::to_string(gpu_id) + " runtime";
    std::string analysis_label;
    analysis_label += "GPU " + std::to_string(gpu_id) + " analysis";

    statistics.Begin_timer(timer_label);

    std::vector<std::pair<uint32_t, uint32_t >> pairs;
    distributor.Get_pairs(1, 0, pairs);

    mutual_info.resize(num_outputs);

    for (unsigned long i = 0; i < pairs.size(); i += pairs_iter) {
        const unsigned long num_pairs = pairs.size() - i < pairs_iter ? pairs.size() - i : pairs_iter;
        for (int j = 0; j < num_pairs; j++) {
            auxIds[j].x = pairs[i + j].first;
            auxIds[j].y = pairs[i + j].second;
        }
        search.mutualInfo(num_pairs, auxIds, &mutual_info.at(0), minMI, minMIPos, numEntriesWithMI);
    }
    cudaDeviceSynchronize();

    statistics.End_timer(timer_label);
    statistics.Add_value(analysis_label, myTotalAnal);
}