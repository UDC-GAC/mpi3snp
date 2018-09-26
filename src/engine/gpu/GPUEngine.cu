/*
 * This file is part of MPI3SNP.
 * Copyright (C) 2018 by Christian Ponte
 *
 * MPI3SNP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MPI3SNP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file gpu/GPUEngine.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief GPUEngine class implementation.
 */

#include "GPUEngine.h"
#include "Dataset.h"
#include "Distributor.h"
#include "EntropySearch.h"
#include <cstring>

GPUEngine::GPUEngine(unsigned int proc_num, unsigned int proc_id,
                     std::vector<std::pair<unsigned int, unsigned int>> gpu_map, bool use_mi, Statistics &statistics) :
        proc_num(proc_num),
        proc_id(proc_id),
        use_mi(use_mi),
        statistics(statistics) {
    statistics.Begin_timer("CUDA initialization time");
    cudaFree(nullptr);
    statistics.End_timer("CUDA initialization time");

    int avail_gpus = 0;
    if (cudaSuccess != cudaGetDeviceCount(&avail_gpus))
        throw CUDAError();
    if (avail_gpus == 0) {
        throw CUDAError("Could not find any CUDA-enabled GPU");
    }

    auto pos = std::find_if(gpu_map.begin(), gpu_map.end(),
                            [&proc_id](std::pair<unsigned int, unsigned int> item) { return item.first == proc_id; });
    gpu_id = pos == gpu_map.end() ? proc_id % avail_gpus : pos->second;

    cudaDeviceProp gpu_prop;
    if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gpu_id))
        throw CUDAError();
    if (gpu_prop.major < 2 || !gpu_prop.canMapHostMemory) {
        throw CUDAError("GPU " + std::to_string(gpu_id) + " does not meet compute capabilities\n" +
                        "Name: " + gpu_prop.name + "\n" + "Compute capability: " +
                        std::to_string(gpu_prop.major) + "." + std::to_string(gpu_prop.minor));
    }
    if (cudaSuccess != cudaSetDevice(gpu_id))
        throw CUDAError();
}

void GPUEngine::run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs) {
    statistics.Begin_timer("SNPs read time");
    Dataset *dataset;
    try {
        dataset = new Dataset(tped, tfam, Dataset::Transposed);
    } catch (const Dataset::ReadError &error) {
        throw Engine::Error(error.what());
    }
    statistics.End_timer("SNPs read time");

    statistics.Addi("SNP count", dataset->Get_SNP_count());
    statistics.Addi("Number of cases", dataset->Get_case_count());
    statistics.Addi("Number of controls", dataset->Get_ctrl_count());

    Distributor<uint32_t, uint2> distributor(dataset->Get_SNP_count(), proc_num);

    EntropySearch search(use_mi, dataset->Get_SNP_count(), dataset->Get_case_count(), dataset->Get_ctrl_count(),
                         dataset->Get_cases(), dataset->Get_ctrls());

    std::vector<uint2> pairs;
    distributor.get_pairs([](uint32_t x, uint32_t y) {
        uint2 p {x, y};
        return p;
    }, proc_id, pairs);

    long myTotalAnal = 0;
    const unsigned int num_snps = dataset->Get_SNP_count();
    for (auto p : pairs) {
        myTotalAnal += num_snps - p.y - 1;
    }
    statistics.Addl("GPU " + std::to_string(gpu_id) + " computations", myTotalAnal);

    std::string timer_label;
    timer_label += "GPU " + std::to_string(gpu_id) + " runtime";
    statistics.Begin_timer(timer_label);

    mutual_info.resize(num_outputs);
    search.mutualInfo(pairs, num_outputs, &mutual_info.at(0));
    cudaDeviceSynchronize();

    delete dataset;

    statistics.End_timer(timer_label);
}