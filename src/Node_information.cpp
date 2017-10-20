//
// Created by christian on 18/10/17.
//

#include "Node_information.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <mpi.h>
#include "Definitions.h"

#ifdef MPI3SNP_USE_GPU

#include "gpu/Gpu_node_information.h"

#else

#include "Cpu_node_information.h"

#endif

Node_information *Node_information::Builder::get_information() {
#ifdef MPI3SNP_USE_GPU
    return new Gpu_node_information();
#else
    return new Cpu_node_information();
#endif
}

Node_information *Node_information::build_from_byteblock(const void *ptr) {
#ifdef MPI3SNP_USE_GPU
    return new Gpu_node_information(ptr);
#else
    return new Cpu_node_information(ptr);
#endif
}

Node_information::Node_information() {
    std::ifstream ifs(hardware_id_file);
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    std::hash<std::string> hash_function;
    hardware_id = std::to_string(hash_function(content));
}

std::vector<Node_information *> Node_information::gather(int process) {
    int rank, count;

    MPI_Comm_size(MPI_COMM_WORLD, &count);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Node_information *local = Node_information::Builder::get_information();
    std::vector<Node_information *> output;

    if (rank == process) {
        output.push_back(local);

        MPI_Status status;
        char *block;
        int block_size;
        Node_information *temp;
        for (int i = 0; i < count; i++) {
            if (i != process) {
                MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_BYTE, &block_size);
                block = new char[block_size];
                MPI_Recv(block, block_size, MPI_BYTE, i, 0, MPI_COMM_WORLD, nullptr);
                temp = Node_information::build_from_byteblock(block);
                delete[] block;
                auto pos = std::find_if(output.begin(), output.end(),
                                        [&temp](Node_information *it) { return it->hardware_id == temp->hardware_id; });
                if (pos == output.end()) {
                    output.push_back(temp);
                } else {
                    (*pos)->add_processes(temp->processes());
                    delete temp;
                }
            }
        }
    } else {
        char *block;
        int size = local->to_byteblock((void **) &block);
        MPI_Send(block, size, MPI_BYTE, process, 0, MPI_COMM_WORLD);
        delete[] block;
        delete local;
    }
    return output;
}
