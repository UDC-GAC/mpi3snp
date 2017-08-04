/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include "Search.h"
#include "Definitions.h"

#ifdef MPI3SNP_USE_GPU

#include "gpu/GPUEngine.h"
#include "gpu/CUDAError.h"

#else

#include "cpu/CPUEngine.h"

#endif

#include <fstream>
#include <thread>
#include "Statistics.h"
#include "IOMpi.h"
#include "Dataset.h"

Search::Builder::Builder(std::string tped_file, std::string tfam_file, std::string out_file) {
    search_obj = new Search();
    search_obj->tped_file = tped_file;
    search_obj->tfam_file = tfam_file;
    search_obj->out_file = out_file;
    search_obj->num_outputs = 10;
    search_obj->use_mi = true;
}

Search::Builder &Search::Builder::Set_num_outputs(unsigned int num_outputs) {
    search_obj->num_outputs = num_outputs;
    return *this;
}

Search::Builder &Search::Builder::Set_gpu_ids(std::vector<unsigned int> gpu_ids) {
    search_obj->gpu_ids = gpu_ids;
    return *this;
}

Search::Builder &Search::Builder::Set_cpu_threads(unsigned int threads) {
    search_obj->cpu_threads = threads;
    return *this;
}

Search::Builder &Search::Builder::Set_use_mi(bool use_mi) {
    search_obj->use_mi = use_mi;
    return *this;
}

Search *Search::Builder::Create_object() {
    return search_obj;
}

Search::Search() :
        cpu_threads(std::thread::hardware_concurrency()) {}

void Search::execute() {
    int proc_id, num_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    std::vector<MutualInfo> mutual_info, result;
    Statistics statistics;

#ifdef MPI3SNP_USE_GPU
    try {
        GPUEngine gpu_engine((unsigned int) num_proc, (unsigned int) proc_id, use_mi);
        gpu_engine.run(tped_file, tfam_file, mutual_info, num_outputs, statistics);
    } catch (const Dataset::ReadError &e) {
        IOMpi::Instance().Mfprintf<IOMpi::E>(std::cerr, (std::string(e.what()) + "\n").c_str());
        MPI_Abort(MPI_COMM_WORLD, 1);
    } catch (const CUDAError &e) {
        IOMpi::Instance().Mfprintf<IOMpi::E>(std::cerr, (std::string(e.what()) + "\n").c_str());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#else
    try {
        CPUEngine cpu_engine(num_proc, proc_id, cpu_threads, use_mi);
        cpu_engine.execute(tped_file, tfam_file, mutual_info, num_outputs, statistics);
    } catch (const Dataset::ReadError &e) {
        IOMpi::Instance().Mfprintf<IOMpi::E>(std::cerr, (std::string(e.what()) + "\n").c_str());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif

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
        for (int i = 0; i < num_outputs; i++) {
            of << result[i].To_string() << '\n';
        }
        of.close();
    }

    IOMpi::Instance().Cprintf<IOMpi::D>("3-SNP analysis finalized\n");
    // Print runtime statistics to stdout
    IOMpi::Instance().Cprintf<IOMpi::B>(statistics.To_string().c_str());
}