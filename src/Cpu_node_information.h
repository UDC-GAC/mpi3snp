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
 * @file Cpu_node_information.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Declaration of the abstract class Node_information, specific for CPU architectures. Provides the host name,
 * cpu name, memory size, specific mpi library implementation and the list of processes hosted in that node.
 */

#ifndef MPI3SNP_CPU_NODE_INFORMATION_H
#define MPI3SNP_CPU_NODE_INFORMATION_H

#include "Node_information.h"

class Cpu_node_information : public Node_information {
public:
    Cpu_node_information();

    explicit Cpu_node_information(const void *ptr);

    std::string hostname() const override;

    std::string cpu() const override;

    long memory() const override;

    std::string mpi_library_version() const override;

    std::vector<int> processes() const override;

    std::vector<std::string> gpus() const override;

    std::string to_string() const override;

protected:
    size_t to_byteblock(void **ptr) const override;

    void add_processes(std::vector<int> processes) override;

private:
    static constexpr const char *cpu_file = "/proc/cpuinfo";

    static constexpr const char *memory_file = "/proc/meminfo";

    static std::string get_cpu_info();

    static long get_physical_memory();

    std::string hostname_str;
    std::string cpu_str;
    long memory_size;
    std::string mpi_library;
    std::vector<int> process_list;
};


#endif //MPI3SNP_CPU_NODE_INFORMATION_H
