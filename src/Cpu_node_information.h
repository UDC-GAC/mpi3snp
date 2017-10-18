//
// Created by christian on 18/10/17.
//

#ifndef MPI3SNP_CPU_NODE_INFORMATION_H
#define MPI3SNP_CPU_NODE_INFORMATION_H


#include "Node_information.h"

class Cpu_node_information : public Node_information {
public:
    Cpu_node_information();

    explicit Cpu_node_information(const void *ptr);

    std::string mpi_library_version() override;

    std::vector<int> processes() override;

protected:
    size_t to_byteblock(void **ptr) override;

private:
    std::string mpi_library;
    std::vector<int> process_list;
};


#endif //MPI3SNP_CPU_NODE_INFORMATION_H
