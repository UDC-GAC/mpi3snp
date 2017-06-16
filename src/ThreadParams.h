/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "MutualInfo.h"
#include "Dataset.h"
#include "Distributor.h"

struct ThreadParams {
    ThreadParams(unsigned int gpu_id, unsigned int num_outputs, Dataset &dataset, Distributor &dist, bool use_MI,
                 Statistics &stats) :
            gpu_id(gpu_id),
            num_outputs(num_outputs),
            dataset(dataset),
            distributor(dist),
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
    Distributor &distributor;
    bool mi;
    // Return values
    Statistics &statistics;
    MutualInfo *mutual_info = nullptr;
};


#endif /* THREADPARAMS_H_ */
