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
    mpi_library.copy(mpi_v, (unsigned long) len-1);
}

Cpu_node_information::Cpu_node_information(const void *ptr) : Node_information() {
    size_t hid_length = *((size_t *) &ptr[0]);
    hardware_id.copy((char *) &ptr[sizeof(size_t)], hid_length);
    size_t mpi_library_len = *((size_t *) &ptr[sizeof(size_t) + hid_length]);
    mpi_library.copy((char *) &ptr[2 * sizeof(size_t) + hid_length], mpi_library_len);
    size_t process_list_size = *((size_t *) &ptr[2 * sizeof(size_t) + hid_length + mpi_library_len]);
    process_list.reserve(process_list_size);
    memcpy(&process_list[0], &ptr[3 * sizeof(size_t) + hid_length + mpi_library_len], process_list_size * sizeof(int));
}

std::string Cpu_node_information::mpi_library_version() {
    return mpi_library;
}

std::vector<int> Cpu_node_information::processes() {
    return std::vector<int>();
}

size_t Cpu_node_information::to_byteblock(void **ptr) {
    size_t hid_size = sizeof(size_t) + hardware_id.length() * sizeof(char);
    size_t mpi_library_size = sizeof(size_t) + mpi_library.length() * sizeof(char);
    size_t process_vector_size = sizeof(size_t) + process_list.size() * sizeof(int);

    *ptr = new char[hid_size + mpi_library_size + process_vector_size];
    *((size_t *) &(*ptr)[sizeof(size_t)]) = hardware_id.length();
    memcpy(&(*ptr)[hid_size], &hardware_id[0], hardware_id.length());
    *((size_t *) &(*ptr)[hid_size]) = mpi_library.length();
    memcpy(&(*ptr)[hid_size + sizeof(size_t)], &mpi_library[0], mpi_library.length());
    *((size_t *) &(*ptr)[hid_size + mpi_library_size]) = process_list.size();
    memcpy(&(*ptr)[hid_size + mpi_library_size + sizeof(size_t)],
           &process_list[0], process_list.size() * sizeof(int));

    return hid_size + mpi_library_size + process_vector_size;
}
