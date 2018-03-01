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
 * @file Node_information.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Abstract class definition. Node_information provides debug information about the nodes where MPI is running,
 * and the process allocation.
 */

#ifndef MPI3SNP_NODE_INFORMATION_H
#define MPI3SNP_NODE_INFORMATION_H

#include <vector>
#include <string>

class Node_information {
public:
    class Builder {
    public:
        static Node_information *get_information();
    };

    Node_information();

    virtual std::string hostname() const =0;

    virtual std::string cpu() const =0;

    virtual long memory() const =0;

    virtual std::string mpi_library_version() const =0;

    virtual std::vector<int> processes() const =0;

    virtual std::vector<std::string> gpus() const =0;

    virtual std::string to_string() const =0;

    static std::vector<Node_information *> gather(int process);

protected:
    virtual size_t to_byteblock(void **ptr) const =0;

    virtual void add_processes(std::vector<int> processes) =0;

    static Node_information *build_from_byteblock(const void *ptr);

    std::string hardware_id;

private:
    static constexpr const char *hardware_id_file = "/proc/net/netlink";
};


#endif //MPI3SNP_NODE_INFORMATION_H
