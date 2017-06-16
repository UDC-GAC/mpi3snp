/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include "GPUSearchMI.h"
#include "GPUEngine.h"
#include "IOMpi.h"
#include "Dataset.h"

GPUSearchMI::Builder::Builder(std::string tped_file, std::string tfam_file, std::string out_file) {
    search_obj = new GPUSearchMI();
    search_obj->tped_file = tped_file;
    search_obj->tfam_file = tfam_file;
    search_obj->out_file = out_file;
    search_obj->num_outputs = 10;
    search_obj->use_mi = true;
}

GPUSearchMI::Builder &GPUSearchMI::Builder::Set_num_outputs(unsigned int num_outputs) {
    search_obj->num_outputs = num_outputs;
    return *this;
}

GPUSearchMI::Builder &GPUSearchMI::Builder::Set_gpu_ids(std::vector<unsigned int> gpu_ids) {
    search_obj->gpu_ids = gpu_ids;
    return *this;
}

GPUSearchMI::Builder &GPUSearchMI::Builder::Set_use_mi(bool use_mi) {
    search_obj->use_mi = use_mi;
    return *this;
}

GPUSearchMI *GPUSearchMI::Builder::Create_object() {
    return search_obj;
}

GPUSearchMI::GPUSearchMI() {
    MPI_Type_contiguous(sizeof(MutualInfo), MPI_CHAR, &MPI_MUTUAL_INFO);
    MPI_Type_commit(&MPI_MUTUAL_INFO);
}

GPUSearchMI::~GPUSearchMI() {
    MPI_Type_free(&MPI_MUTUAL_INFO);
}

void GPUSearchMI::execute() {
    int proc_id, num_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    std::vector<MutualInfo> mutual_info;
    Statistics statistics;

    try {
        GPUEngine gpu_engine((unsigned int) num_proc, (unsigned int) proc_id, gpu_ids, use_mi);
        gpu_engine.run(tped_file, tfam_file, mutual_info, num_outputs, statistics);

        if (proc_id == 0) {
            // Gather all the results
            std::vector<MutualInfo> buff(num_outputs);
            for (int rank = 1; rank < num_proc; rank++) {
                MPI_Recv(&buff[0], num_outputs, MPI_MUTUAL_INFO, rank, 123, MPI_COMM_WORLD, NULL);
                mutual_info.insert(mutual_info.end(), buff.begin(), buff.end());
            }

            // Sort the result
            std::sort(&mutual_info[0], &mutual_info[num_proc * num_outputs]);

            // Write results to the output file
            FILE *out = fopen(out_file.c_str(), "wb");
            if (out == NULL) {
                IOMpi::Instance().Mprintf("Out file: file %s could not be opened\n", out_file.c_str());
                return;
            }
            std::vector<MutualInfo>::reverse_iterator it;
            int count;
            for (it = mutual_info.rbegin(), count = 0; count < num_outputs; it++, count++) {
                fprintf(out, "%u %u %u %f\n", it->_id1, it->_id2, it->_id3, it->_mutualInfoValue);
            }
            fclose(out);
        } else {
            // Send results to master
            MPI_Send(&mutual_info[0], num_outputs, MPI_MUTUAL_INFO, 0, 123, MPI_COMM_WORLD);
        }

        IOMpi::Instance().Cprintf("3-SNP analysis finalized\n");

        // Print runtime statistics to stdout
        auto timers = statistics.Get_all_timers();
        std::string output("Statistics\n");
        for (auto it = timers.begin(); it < timers.end(); it++) {
            output += "\t" + it->first + ": " + std::to_string(it->second) + " seconds\n";
        }
        auto values = statistics.Get_all_values();
        for (auto it = values.begin(); it < values.end(); it++) {
            output += "\t" + it->first + ": " + std::to_string(it->second) + "\n";
        }
        IOMpi::Instance().Cprintf(output.c_str());
    } catch (Dataset::ReadError &e) {
        IOMpi::Instance().Mprintf((std::string(e.what()) + "\n").c_str());
    }
}