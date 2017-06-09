/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "Macros.h"
#include "GPUSNPDistributorStatic.h"
#include "GPUEngine.h"
#include "MutualInfo.h"

struct ThreadParams {
    ThreadParams(int tid, uint16_t numOutputs, GPUSNPDistributor *distributor, int gpu_id, bool isMI,
                 Statistics &stats) : statistics(stats) {
        _tid = tid;
        _numOutputs = numOutputs;
        _distributor = distributor;
        _gpu = gpu_id;
        _isMI = isMI;

        _mutualInfo = new MutualInfo[numOutputs];
    }

    ~ThreadParams() {
        delete[] _mutualInfo;
    }

    // Parameters
    int _tid;
    uint16_t _numOutputs;
    GPUSNPDistributor *_distributor;
    int _gpu;
    bool _isMI;

    // Return values
    Statistics &statistics;
    MutualInfo *_mutualInfo;
};


#endif /* THREADPARAMS_H_ */
