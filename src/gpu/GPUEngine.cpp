/*
 * GPUEngine.cu
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#include "GPUEngine.h"
#include "../Dataset.h"
#include "../Distributor.h"
#include "../IOMpi.h"
#include "EntropySearch.h"
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

    cudaDeviceProp gpu_prop;

    IOMpi::Instance().mprint<IOMpi::D>("Available GPUs:\n");
    for (int gid = 0; gid < avail_gpus; gid++) {
        if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gid))
            throw CUDAError(cudaGetLastError());
        IOMpi::Instance().mprint<IOMpi::D>("GPU " + std::to_string(gid) + ": " + gpu_prop.name + "\n");
    }

    gpu_id = proc_id % avail_gpus;
    if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gpu_id))
        throw CUDAError(cudaGetLastError());
    if (gpu_prop.major < 2 || !gpu_prop.canMapHostMemory) {
        throw CUDAError("GPU " + std::to_string(gpu_id) + " does not meet compute capabilities\n" +
                        "Name: " + gpu_prop.name + "\n" + "Compute capability: " +
                        std::to_string(gpu_prop.major) + "." + std::to_string(gpu_prop.minor));
    }
    if (cudaSuccess != cudaSetDevice(gpu_id))
        throw CUDAError(cudaGetLastError());

    IOMpi::Instance().print<IOMpi::D>("Process " + std::to_string(proc_id) +
                                      " using GPU " + std::to_string(gpu_id) + "\n");
}

void GPUEngine::run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs,
                    Statistics &statistics) {
    statistics.Begin_timer("SNPs read time");
    Dataset dataset(tped, tfam, Dataset::Transposed);
    statistics.End_timer("SNPs read time");
    statistics.Add_value("SNP count", dataset.Get_SNP_count());
    statistics.Add_value("Number of cases", dataset.Get_case_count());
    statistics.Add_value("Number of controls", dataset.Get_ctrl_count());

    Distributor distributor(proc_num, proc_id, dataset.Get_SNP_count());

    EntropySearch search(use_mi, dataset.Get_SNP_count(), dataset.Get_case_count(), dataset.Get_ctrl_count(),
                         dataset.Get_cases(), dataset.Get_ctrls());

    std::vector<std::pair<uint32_t, uint32_t >> pairs;
    distributor.Get_pairs(1, 0, pairs);

    int myTotalAnal = 0;
    const unsigned int num_snps = dataset.Get_SNP_count();
    for (auto p : pairs) {
        myTotalAnal += num_snps - p.second - 1;
    }
    statistics.Add_value("GPU " + std::to_string(gpu_id) + " computations", myTotalAnal);

    std::string timer_label;
    timer_label += "GPU " + std::to_string(gpu_id) + " runtime";
    statistics.Begin_timer(timer_label);

    mutual_info.resize(num_outputs);
    search.mutualInfo(pairs, num_outputs, &mutual_info.at(0));
    cudaDeviceSynchronize();

    statistics.End_timer(timer_label);
}