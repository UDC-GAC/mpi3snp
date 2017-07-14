/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "../MutualInfo.h"
#include "../Dataset.h"
#include "../Distributor.h"

struct ThreadParams {
    ThreadParams(unsigned int gpu_id, unsigned int num_outputs, Dataset &dataset, bool use_MI, Statistics &stats) :
            gpu_id(gpu_id),
            num_outputs(num_outputs),
            dataset(dataset),
            mi(use_MI),
            statistics(stats) {
        mutual_info = new MutualInfo[num_outputs];
    }

    ~ThreadParams() {
        delete[] mutual_info;
    }

    // Parameters
    unsigned int gpu_id;
    unsigned int num_outputs;
    Dataset &dataset;
    std::vector<std::pair<uint32_t, uint32_t >> pairs;
    bool mi;
    // Return values
    Statistics &statistics;
    MutualInfo *mutual_info = nullptr;
};


#endif /* THREADPARAMS_H_ */
