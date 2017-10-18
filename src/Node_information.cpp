//
// Created by christian on 18/10/17.
//

#include "Node_information.h"
#include <array>
#include <memory>
#include <mpi.h>
#include "Definitions.h"
#include "Cpu_node_information.h"

Node_information* Node_information::Builder::get_information() {
#ifdef MPI3SNP_USE_GPU
    return new Gpu_node_information();
#else
    return new Cpu_node_information();
#endif
}

Node_information* Node_information::Builder::build_from_byteblock(const void *ptr) {
#ifdef MPI3SNP_USE_GPU
    return new Gpu_node_information();
#else
    return new Cpu_node_information(ptr);
#endif
}

Node_information::Node_information() {
    hardware_id.clear();
    std::array<char, 128> buffer;
    std::shared_ptr<FILE> pipe(popen(hid_bash_command, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            hardware_id.append(buffer.data());
    }
}

std::vector<Node_information> Node_information::gather(int process) {
    std::vector<Node_information> output;
    int proc_num;

    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    output.reserve((unsigned long) proc_num);
    Node_information *local = Node_information::Builder::get_information();
    void *byteblock;
    size_t blocksize = local->to_byteblock(&byteblock);
    delete local;



    return std::vector<Node_information>();
}
