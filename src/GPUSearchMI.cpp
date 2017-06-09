/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include "GPUSearchMI.h"
#include "GPUEngine.h"
#include "IOMpi.h"

GPUSearchMI::GPUSearchMI(Options *options) : Search(options) {
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

    GPUEngine gpu_engine(_options);
    gpu_engine.run(mutual_info, statistics);

    if (proc_id == 0) {
        // Gather all the results
        std::vector<MutualInfo> buff(_options->getNumOutputs());
        for (int rank = 1; rank < num_proc; rank++) {
            MPI_Recv(&buff[0], _options->getNumOutputs(), MPI_MUTUAL_INFO, rank, 123, MPI_COMM_WORLD, NULL);
            mutual_info.insert(mutual_info.end(), buff.begin(), buff.end());
        }

        // Sort the result
        std::sort(&mutual_info[0], &mutual_info[num_proc * _options->getNumOutputs()]);

        // Write results to the output file
        MyFilePt out = myfopen(_options->getOutFileName(), "wb");
        if (out == NULL) {
            IOMpi::Instance().Mprintf("Out file: file %s could not be opened\n", _options->getOutFileName());
            return;
        }
        LineReader lr;
        std::vector<MutualInfo>::reverse_iterator it;
        int count;
        for (it = mutual_info.rbegin(), count = 0; count < _options->getNumOutputs(); it++, count++) {
            lr.printMI(out, it->_id1, it->_id2, it->_id3, it->_mutualInfoValue);
        }
        myfclose(out);
    } else {
        // Send results to master
        MPI_Send(&mutual_info[0], _options->getNumOutputs(), MPI_MUTUAL_INFO, 0, 123, MPI_COMM_WORLD);
    }

    IOMpi::Instance().Cprintf("3-SNP analysis finalized\n");

    // Print runtime statistics to stdout
    auto timers = statistics.Get_all_timers();
    std::string output("Statistics\n");
    for (auto it = timers.begin(); it < timers.end(); it++) {
        output += "\t" + it->first + ": " + std::to_string(it->second) + " seconds\n";
    }
    auto values = statistics.Get_all_values();
    for (auto it = values.begin(); it < values.end(); it++){
        output += "\t" + it->first + ": " + std::to_string(it->second) + "\n";
    }
    IOMpi::Instance().Cprintf(output.c_str());

}