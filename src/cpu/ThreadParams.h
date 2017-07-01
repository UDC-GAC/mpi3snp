/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "Macros.h"
#include "../MutualInfo.h"
#include "../Distributor.h"
#include "../Dataset.h"

struct ThreadParams {
    ThreadParams(int tid, Distributor &distributor, Dataset &dataset, uint16_t numOutputs) :
            tid(tid),
            distributor(distributor),
            dataset(dataset),
            numOutputs(numOutputs),
            mutualInfo(new MutualInfo[numOutputs]) {
        _numAnalyzed = 0;
        _runtime = 0;
    }

    ~ThreadParams() {
        delete[] mutualInfo;
    }

    const int tid;
    Distributor &distributor;
    Dataset &dataset;
    const uint16_t numOutputs;
    MutualInfo *mutualInfo;
    // Values for the return information
    uint64_t _numAnalyzed;
    double _runtime;
};


#endif /* THREADPARAMS_H_ */
