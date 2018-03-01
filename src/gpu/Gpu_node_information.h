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
 * @file gpu/Gpu_node_information.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief GPU implementation of the abstract class Node_information, extending the CPU implementation provided in
 * Cpu_node_information. This class adds information about the GPUs present in the node.
 */

#ifndef MPI3SNP_GPU_NODE_INFORMATION_H
#define MPI3SNP_GPU_NODE_INFORMATION_H

#include "../Cpu_node_information.h"

class Gpu_node_information : public Cpu_node_information {
public:
    Gpu_node_information();

    explicit Gpu_node_information(const void *ptr);

    std::vector<std::string> gpus() const override;

    std::string to_string() const override;

protected:
    size_t to_byteblock(void **ptr) const override;

private:
    std::vector<std::string> gpu_list;
};

#endif //MPI3SNP_GPU_NODE_INFORMATION_H
