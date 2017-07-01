/*
 * Engine.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#include "CPUEngine.h"
#include "ThreadParams.h"
#include "EntropySearch.h"
#include "../ThreadError.h"
#include "../Dataset.h"
#include "../Distributor.h"
#include "../IOMpi.h"
#include <cfloat>
#include <vector>

CPUEngine::CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi) {
    this->num_proc = num_proc;
    this->proc_id = proc_id;
    this->num_threads = num_threads;
    this->use_mi = use_mi;
}

CPUEngine::~CPUEngine() {

}

void CPUEngine::execute(std::string tped_file, std::string tfam_file, std::vector<MutualInfo> &mutual_info,
                        uint16_t num_outputs, Statistics &statistics) {
    double stime = MPI_Wtime();
    double etime;

    Dataset dataset(tped_file, tfam_file, Dataset::Regular);

    Distributor distributor(num_proc, proc_id, dataset.Get_SNP_count(), num_threads > 0);

    etime = MPI_Wtime();

    IOMpi::Instance().Cprintf<IOMpi::D>("Loaded %ld SNPs (%ld/%ld cases/controls) in %.2f seconds\n",
                                        dataset.Get_SNP_count(), dataset.Get_case_count(), dataset.Get_ctrl_count(),
                                        etime - stime);

    std::vector<pthread_t> threadIDs(num_threads, 0);

    // Computation of the single-SNP entropy
    std::vector<ThreadParams *> params(num_threads);
    for (int tid = 0; tid < num_threads; tid++) {
        params[tid] = new ThreadParams(tid, distributor, dataset, num_outputs);

        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], NULL, threadMI, params[tid]) != 0) {
            for (int i = 0; i < tid; i++)
                pthread_cancel(threadIDs[i]);
            char message[80];
            sprintf(message, "error while creating the thread %i in process %i", tid, proc_id);
            throw ThreadError(message);
        }
    }

    MutualInfo *auxMutualInfo = new MutualInfo[num_threads * num_outputs];

    // Wait for the completion of all threads
    for (int tid = 0; tid < num_threads; tid++) {
        pthread_join(threadIDs[tid], NULL);
        memcpy(&auxMutualInfo[tid * num_outputs], params[tid]->mutualInfo, num_outputs * sizeof(MutualInfo));
        delete params[tid];
    }

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo, auxMutualInfo + num_outputs * num_threads);
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], auxMutualInfo + num_outputs * (num_threads - 1), sizeof(MutualInfo) * num_outputs);
}

void *CPUEngine::threadMI(void *arg) {
    // Enable thread cancellation
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, nullptr);
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, nullptr);

    ThreadParams *params = (ThreadParams *) arg;

    EntropySearch search(params->dataset.Get_SNP_count(), params->dataset.Get_case_count(), params->dataset.Get_cases(),
                         params->dataset.Get_ctrl_count(), params->dataset.Get_ctrls());

    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    pairs.reserve(Distributor::DEFAULT_PAIRS_BLOCK);

    // Variables to work with the outputs
    MutualInfo *mutualInfo = new MutualInfo[params->numOutputs];
    // The minimum value in the array
    float minMI = FLT_MAX;
    // The position of the minimum value
    uint16_t minMIPos = 0;
    // Number of entries of the array full
    uint16_t numEntriesWithMI = 0;

    uint64_t myTotalAnal = 0;
    double stime = MPI_Wtime();

    while (params->distributor.Get_pairs(pairs, myTotalAnal) > 0) {
        search.mutualInfo(pairs, mutualInfo, params->numOutputs, minMI, minMIPos, numEntriesWithMI);
        pairs.clear();
    }

    params->_numAnalyzed = myTotalAnal;
    params->_runtime = MPI_Wtime() - stime;
    memcpy(params->mutualInfo, mutualInfo, params->numOutputs * sizeof(MutualInfo));

    delete[] mutualInfo;
    return NULL;
}

