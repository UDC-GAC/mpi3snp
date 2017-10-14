/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include "Search.h"
#include <fstream>
#include "IOMpi.h"
#include "Dataset.h"
#include "Definitions.h"

#ifdef MPI3SNP_USE_GPU

#include "gpu/GPUEngine.h"
#include "gpu/CUDAError.h"

#else

#include "cpu/CPUEngine.h"

#endif

Search *Search::Builder::build_from_args(Arg_parser::Arguments arguments) {
    auto *search = new Search();
    MPI_Comm_rank(MPI_COMM_WORLD, &search->proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &search->num_proc);
    search->tped_file = arguments.tped;
    search->tfam_file = arguments.tfam;
    search->out_file = arguments.output;
    search->num_outputs = arguments.output_num;

#ifdef MPI3SNP_USE_GPU
    search->engine = new GPUEngine(search->num_proc, search->proc_id, arguments.gpu_map, arguments.use_mi);
#else
    search->engine = new CPUEngine(search->num_proc, search->proc_id, arguments.cpu_threads, arguments.use_mi);
#endif
    return search;
}

void Search::execute() {
    std::vector<MutualInfo> mutual_info, result;
    Statistics statistics;

    try {
        engine->run(tped_file, tfam_file, mutual_info, num_outputs, statistics);
    } catch (const Dataset::ReadError &e) {
        IOMpi::Instance().smprint<IOMpi::E>(std::cerr, std::string(e.what()) + "\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
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
    // Print runtime statistics to stdout
    IOMpi::Instance().cprint<IOMpi::B>(statistics.To_string());
}