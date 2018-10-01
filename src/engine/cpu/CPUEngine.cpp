/*
 * MPI3SNP ~ https://github.com/chponte/mpi3snp
 *
 * Copyright 2018 Christian Ponte
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file cpu/CPUEngine.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief CPUEngine class members implementation.
 */

#include "CPUEngine.h"
#include "ThreadParams.h"
#include "EntropySearch.h"
#include "ThreadError.h"
#include <cfloat>

CPUEngine::CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi, Statistics &statistics) :
        num_proc(num_proc),
        proc_id(proc_id),
        num_threads(num_threads),
        use_mi(use_mi),
        statistics(statistics) {}

void CPUEngine::run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs) {
    statistics.Begin_timer("SNPs read time");
    Dataset *dataset;
    try {
        dataset = new Dataset(tped, tfam, Dataset::Regular);
    } catch (const Dataset::Read_error &error) {
        throw Engine::Error(error.what());
    }
    statistics.End_timer("SNPs read time");

    Distributor<uint32_t, std::pair<uint32_t, uint32_t>> distributor(dataset->get_SNP_count(), num_proc * num_threads);
    statistics.Addi("SNP count", dataset->get_SNP_count());
    statistics.Addi("Number of cases", dataset->get_case_count());
    statistics.Addi("Number of controls", dataset->get_ctrl_count());

    std::vector<pthread_t> threadIDs(num_threads, 0);

    // Computation of the single-SNP entropy
    std::vector<ThreadParams *> params(num_threads);
    for (int tid = 0; tid < num_threads; tid++) {
        params[tid] = new ThreadParams(tid, *dataset, num_outputs, statistics);
        distributor.get_pairs([](uint32_t first, uint32_t second) { return std::make_pair(first, second); },
                              proc_id * num_threads + tid, params[tid]->pairs);

        // Create thread entities that call to the functions below
        if (pthread_create(&threadIDs[tid], nullptr, thread, params[tid]) != 0) {
            for (int i = 0; i < tid; i++)
                pthread_cancel(threadIDs[i]);
            throw ThreadError("Error creating thread number " + std::to_string(tid));
        }
    }

    MutualInfo auxMutualInfo[num_threads * num_outputs];

    // Wait for the completion of all threads
    for (int tid = 0; tid < num_threads; tid++) {
        pthread_join(threadIDs[tid], nullptr);
        memcpy(&auxMutualInfo[tid * num_outputs], params[tid]->mutualInfo, num_outputs * sizeof(MutualInfo));
        delete params[tid];
    }

    delete dataset;

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo, auxMutualInfo + num_outputs * num_threads);
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], auxMutualInfo + num_outputs * (num_threads - 1), sizeof(MutualInfo) * num_outputs);
}

void *CPUEngine::thread(void *arg) {
    // Enable thread cancellation
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, nullptr);
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, nullptr);

    ThreadParams *params = (ThreadParams *) arg;

    EntropySearch search(params->dataset.get_SNP_count(), params->dataset.get_case_count(), params->dataset.get_cases(),
                         params->dataset.get_ctrl_count(), params->dataset.get_ctrls());

    // Variables to work with the outputs
    MutualInfo *mutualInfo = new MutualInfo[params->numOutputs];
    // The minimum value in the array
    float minMI = FLT_MAX;
    // The position of the minimum value
    uint16_t minMIPos = 0;
    // Number of entries of the array full
    uint16_t numEntriesWithMI = 0;

    long myTotalAnal = 0;
    std::string timer_label;
    timer_label += "Thread " + std::to_string(params->tid) + " runtime";
    std::string analysis_label;
    analysis_label += "Thread " + std::to_string(params->tid) + " computations";

    params->statistics.Begin_timer(timer_label);

    myTotalAnal = search.mutualInfo(params->pairs, mutualInfo, params->numOutputs, minMI, minMIPos, numEntriesWithMI);

    params->statistics.End_timer(timer_label);
    params->statistics.Addl(analysis_label, myTotalAnal);

    memcpy(params->mutualInfo, mutualInfo, params->numOutputs * sizeof(MutualInfo));

    delete[] mutualInfo;
    return nullptr;
}

