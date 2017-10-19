//
// Created by christian on 18/10/17.
//

#include "Cpu_node_information.h"
#include <cstring>
#include <mpi.h>

Cpu_node_information::Cpu_node_information() : Node_information() {
    int id, len;
    char mpi_v[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    process_list.push_back(id);
    MPI_Get_library_version(mpi_v, &len);
    mpi_library.append(mpi_v, (unsigned long) len - 1);
}

Cpu_node_information::Cpu_node_information(const void *pointer) : Node_information() {
    auto *ptr = (char *) pointer;
    size_t size, offset = 0;

    // Hardware ID string
    memcpy(&size, ptr, sizeof(size_t));
    offset += sizeof(size_t);
    hardware_id.resize(size);
    memcpy(&hardware_id[0], ptr + offset, size);
    offset += size;

    // MPI library string
    memcpy(&size, ptr + offset, sizeof(size_t));
    offset += sizeof(size_t);
    mpi_library.resize(size);
    memcpy(&mpi_library[0], ptr + offset, size);
    offset += size;

    // Process ID vector
    memcpy(&size, ptr + offset, sizeof(size_t));
    offset += sizeof(size_t);
    process_list.resize(size);
    memcpy(&process_list[0], ptr + offset, size * sizeof(int));
    offset += size * sizeof(int);
}

std::string Cpu_node_information::mpi_library_version() const {
    return mpi_library;
}

std::vector<int> Cpu_node_information::processes() const {
    return process_list;
}

size_t Cpu_node_information::to_byteblock(void **pointer) const {
    size_t hid_size = sizeof(size_t) + hardware_id.length();
    size_t mpi_library_size = sizeof(size_t) + mpi_library.length();
    size_t process_vector_size = sizeof(size_t) + process_list.size() * sizeof(int);
    *pointer = new char[hid_size + mpi_library_size + process_vector_size];
    auto *ptr = (char *) *pointer;
    size_t temp, offset = 0;

    // Hardware ID string
    temp = hardware_id.length();
    memcpy(ptr + offset, &temp, sizeof(size_t));
    offset += sizeof(size_t);
    memcpy(ptr + offset, &hardware_id[0], hardware_id.length());
    offset += hardware_id.length();

    // MPI library string
    temp = mpi_library.length();
    memcpy(ptr + offset, &temp, sizeof(size_t));
    offset += sizeof(size_t);
    memcpy(ptr + offset, &mpi_library[0], mpi_library.length());
    offset += mpi_library.length();

    // Process ID vector
    temp = process_list.size();
    memcpy(ptr + offset, &temp, sizeof(size_t));
    offset += sizeof(size_t);
    memcpy(ptr + offset, &process_list[0], process_list.size() * sizeof(int));
    offset += process_list.size() * sizeof(int);

    return offset;
}

void Cpu_node_information::add_processes(std::vector<int> processes) {
    process_list.insert(process_list.end(), processes.begin(), processes.end());
}
