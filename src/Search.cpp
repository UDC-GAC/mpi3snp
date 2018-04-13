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
 * @file Search.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Search class members implementation.
 */

#include "Search.h"
#include <fstream>
#include "IOMpi.h"
#include "Dataset.h"
#include "Definitions.h"

#if TARGET_ARCH == CPU

#include "cpu/CPUEngine.h"

#else

#include "gpu/GPUEngine.h"
#include "gpu/CUDAError.h"

#endif

Search *Search::Builder::build_from_args(Arg_parser::Arguments arguments, Statistics &statistics) {
    auto *search = new Search();
    MPI_Comm_rank(MPI_COMM_WORLD, &search->proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &search->num_proc);
    search->tped_file = arguments.tped;
    search->tfam_file = arguments.tfam;
    search->out_file = arguments.output;
    search->num_outputs = arguments.output_num;

#if TARGET_ARCH == CPU
    search->engine = new CPUEngine(search->num_proc, search->proc_id, arguments.cpu_threads, arguments.use_mi,
                                   statistics);
#else
    try {
        search->engine = new GPUEngine(search->num_proc, search->proc_id, arguments.gpu_map, arguments.use_mi, statistics);
    } catch (const Engine::Error &e) {
        IOMpi::Instance().smprint<IOMpi::E>(std::cerr, std::string(e.what()) + "\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return nullptr;
    }
#endif
    return search;
}

void Search::execute() {
    std::vector<MutualInfo> mutual_info, result;

    try {
        engine->run(tped_file, tfam_file, mutual_info, num_outputs);
    } catch (const Engine::Error &e) {
        IOMpi::Instance().smprint<IOMpi::E>(std::cerr, std::string(e.what()) + "\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Gather the results on the master process
    if (proc_id == 0) {
        result.resize(num_outputs * num_proc);
    }

    MPI_Gather(&mutual_info.front(), num_outputs * sizeof(MutualInfo), MPI_BYTE,
               &result.front(), num_outputs * sizeof(MutualInfo), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (proc_id == 0) {
        // Sort the result
        std::sort(result.begin(), result.end(), [](MutualInfo a, MutualInfo b) { return b < a; });

        // Write results to the output file
        std::ofstream of(out_file.c_str(), std::ios::out);
        for (unsigned int i = 0; i < num_outputs; i++) {
            of << result[i].To_string() << '\n';
        }
        of.close();
    }

    IOMpi::Instance().cprint<IOMpi::D>("3-SNP analysis finalized\n");
}