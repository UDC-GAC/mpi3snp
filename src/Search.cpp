/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include <fstream>
#include "Search.h"
#include "gpu/GPUEngine.h"
#include "IOMpi.h"
#include "Dataset.h"
#include "gpu/CUDAError.h"

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

Search::Builder &Search::Builder::Set_use_mi(bool use_mi) {
    search_obj->use_mi = use_mi;
    return *this;
}

Search *Search::Builder::Create_object() {
    return search_obj;
}

Search::Search() {
    MPI_Type_contiguous(sizeof(MutualInfo), MPI_CHAR, &MPI_MUTUAL_INFO);
    MPI_Type_commit(&MPI_MUTUAL_INFO);
}

Search::~Search() {
    MPI_Type_free(&MPI_MUTUAL_INFO);
}

void Search::execute() {
    int proc_id, num_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    std::vector<MutualInfo> mutual_info;
    Statistics statistics;

    try {
        GPUEngine gpu_engine((unsigned int) num_proc, (unsigned int) proc_id, gpu_ids, use_mi);
        gpu_engine.run(tped_file, tfam_file, mutual_info, num_outputs, statistics);
    } catch (const Dataset::ReadError &e) {
        IOMpi::Instance().Mfprintf<IOMpi::E>(std::cerr, (std::string(e.what()) + "\n").c_str());
    } catch (const CUDAError &e) {
        IOMpi::Instance().Mfprintf<IOMpi::E>(std::cerr, (std::string(e.what()) + "\n").c_str());
    }

    if (proc_id == 0) {
        // Gather all the results
        std::vector<MutualInfo> buff(num_outputs);
        for (int rank = 1; rank < num_proc; rank++) {
            MPI_Recv(&buff[0], num_outputs, MPI_MUTUAL_INFO, rank, 123, MPI_COMM_WORLD, NULL);
            mutual_info.insert(mutual_info.end(), buff.begin(), buff.end());
        }

        // Sort the result
        std::sort(&mutual_info[0], &mutual_info[num_proc * num_outputs],
                  [](MutualInfo a, MutualInfo b) { return b < a; });

        // Write results to the output file
        std::ofstream of(out_file.c_str(), std::ios::out);
        for (int i=0; i<num_outputs; i++) {
            of << mutual_info[i].To_string() << '\n';
        }
        of.close();
    } else {
        // Send results to master
        MPI_Send(&mutual_info[0], num_outputs, MPI_MUTUAL_INFO, 0, 123, MPI_COMM_WORLD);
    }

    IOMpi::Instance().Cprintf<IOMpi::D>("3-SNP analysis finalized\n");
    // Print runtime statistics to stdout
    IOMpi::Instance().Cprintf<IOMpi::B>(statistics.To_string().c_str());
}