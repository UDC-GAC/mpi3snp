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
 * @file cpu/CPUEngine.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief CPUEngine class declaration, implementing the abstract class Engine.
 */

#ifndef MPI3SNP_CPUENGINE_H
#define MPI3SNP_CPUENGINE_H

#include "Engine.h"

class CPUEngine : public Engine {
public:
    CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi, Statistics &statistics);

    virtual ~CPUEngine() = default;

    void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs) override;

private:
    static void *threadMI(void *arg);

    int num_proc;
    int proc_id;
    int num_threads;
    bool use_mi;
    Statistics &statistics;
};

#endif //MPI3SNP_CPUENGINE_H
