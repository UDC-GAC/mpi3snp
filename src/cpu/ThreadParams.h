/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "../MutualInfo.h"
#include "../Distributor.h"
#include "../Dataset.h"

struct ThreadParams {
    ThreadParams(int tid, Distributor &distributor, Dataset &dataset, uint16_t numOutputs, Statistics &statistics) :
            tid(tid),
            distributor(distributor),
            dataset(dataset),
            numOutputs(numOutputs),
            mutualInfo(new MutualInfo[numOutputs]),
            statistics(statistics) {}

    ~ThreadParams() {
        delete[] mutualInfo;
    }

    const int tid;
    Distributor &distributor;
    Dataset &dataset;
    const uint16_t numOutputs;
    MutualInfo *mutualInfo;
    Statistics &statistics;
};


#endif /* THREADPARAMS_H_ */
