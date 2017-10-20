//
// Created by christian on 20/10/17.
//

#include "Gpu_node_information.h"
#include <algorithm>
#include <cstring>

Gpu_node_information::Gpu_node_information() : Cpu_node_information() {
    gpu_list.push_back("ejemplo de GPU");
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
        memcpy(block + offset, gpu_info, size);
        offset += size;
        g = std::string(gpu_info);
    }
}

std::vector<std::string> Gpu_node_information::gpus() const {
    return gpu_list;
}

size_t Gpu_node_information::to_byteblock(void **ptr) const {
    void *cpu_info;
    size_t cpu_block_size = Cpu_node_information::to_byteblock(&cpu_info);
    // Size = number of strings + (length of string + string) for each string in the vector
    size_t gpu_vector_size = sizeof(size_t);
    std::for_each(gpu_list.begin(), gpu_list.end(),
                  [&gpu_vector_size](auto str) { gpu_vector_size = +sizeof(size_t) + str.length(); });
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
    for (auto g : gpu_list){
        size = g.length();
        memcpy(buffer + offset, &size, sizeof(size_t));
        offset += sizeof(size_t);
        memcpy(buffer + offset, &g[0], size);
        offset += size;
    }

    return offset;
}
