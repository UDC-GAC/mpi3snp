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
 * @file gpu/GPUEngine.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief GPUEngine class declaration.
 */

#ifndef MPI3SNP_GPUENGINE_H
#define MPI3SNP_GPUENGINE_H

#include "../Engine.h"

class GPUEngine : public Engine {
public:
    GPUEngine(unsigned int proc_num, unsigned int proc_id, std::vector<std::pair<unsigned int, unsigned int>> gpu_map,
              bool use_mi, Statistics &statistics);

    void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs) override;

private:
    unsigned int proc_num, proc_id, gpu_id;
    bool use_mi;
    Statistics &statistics;
};

#endif //MPI3SNP_GPUENGINE_H
