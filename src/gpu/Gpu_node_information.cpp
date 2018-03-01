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
 * @file gpu/Gpu_node_information.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Gpu_node_information class members implementation.
 */

#include "Gpu_node_information.h"
#include <algorithm>
#include <cstring>
#include <cuda_runtime.h>
#include <iostream>

Gpu_node_information::Gpu_node_information() : Cpu_node_information() {
    int avail_gpus = 0;
    if (cudaSuccess != cudaGetDeviceCount(&avail_gpus)) {
        return;
    }

    cudaDeviceProp gpu_prop;
    for (int gid = 0; gid < avail_gpus; gid++) {
        if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gid)) {
            continue;
        }
        gpu_list.push_back(std::to_string(gid) + ": " + gpu_prop.name);
    }
}

Gpu_node_information::Gpu_node_information(const void *ptr) : Cpu_node_information(((char *) ptr) + sizeof(size_t)) {
    auto *block = (char *) ptr;
    size_t size, offset = 0;
    memcpy(&offset, block, sizeof(size_t));

    // GPU information vector
    memcpy(&size, block + offset, sizeof(size_t));
    offset += sizeof(size_t);
    gpu_list.resize(size);
    for (auto &g : gpu_list) {
        memcpy(&size, block + offset, sizeof(size_t));
        offset += sizeof(size_t);
        char gpu_info[size];
        memcpy(gpu_info, block + offset, size);
        offset += size;
        g = std::string(gpu_info, size);
    }
}

std::vector<std::string> Gpu_node_information::gpus() const {
    return gpu_list;
}

std::string Gpu_node_information::to_string() const {
    std::string output = Cpu_node_information::to_string();
    output += "GPUs: ";
    for (const auto &g : gpu_list){
        output += "[" + g + "] ";
    }
    output += "\n";
    return output;
}

size_t Gpu_node_information::to_byteblock(void **ptr) const {
    void *cpu_info;
    size_t cpu_block_size = Cpu_node_information::to_byteblock(&cpu_info);
    // Size = number of strings + (length of string + string) for each string in the vector
    size_t gpu_vector_size = sizeof(size_t);
    std::for_each(gpu_list.begin(), gpu_list.end(),
                  [&gpu_vector_size](auto str) { gpu_vector_size += sizeof(size_t) + str.length(); });
    // Memory allocation
    *ptr = new char[sizeof(size_t) + cpu_block_size + gpu_vector_size];
    auto *buffer = (char *) *ptr;
    size_t size, offset = 0;

    // Length of CPU information byte block
    size = sizeof(size_t) + cpu_block_size;
    memcpy(buffer + offset, &size, sizeof(size_t));
    offset += sizeof(size_t);

    // CPU information byte block
    memcpy(buffer + offset, cpu_info, cpu_block_size);
    offset += cpu_block_size;

    // GPU information vector
    size = gpu_list.size();
    memcpy(buffer + offset, &size, sizeof(size_t));
    offset += sizeof(size_t);
    for (const auto &g : gpu_list) {
        size = g.length();
        memcpy(buffer + offset, &size, sizeof(size_t));
        offset += sizeof(size_t);
        memcpy(buffer + offset, &g[0], size);
        offset += size;
    }
    return offset;
}
