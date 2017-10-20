//
// Created by christian on 20/10/17.
//

#ifndef MPI3SNP_GPU_NODE_INFORMATION_H
#define MPI3SNP_GPU_NODE_INFORMATION_H

#include "../Cpu_node_information.h"

class Gpu_node_information : public Cpu_node_information {
public:
    Gpu_node_information();

    explicit Gpu_node_information(const void *ptr);

    std::vector<std::string> gpus() const override;

protected:
    size_t to_byteblock(void **ptr) const override;

private:
    std::vector<std::string> gpu_list;
};

#endif //MPI3SNP_GPU_NODE_INFORMATION_H
