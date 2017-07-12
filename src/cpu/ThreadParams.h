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
    ThreadParams(int tid, Dataset &dataset, uint16_t numOutputs, Statistics &statistics) :
            tid(tid),
            dataset(dataset),
            numOutputs(numOutputs),
            mutualInfo(new MutualInfo[numOutputs]),
            statistics(statistics) {}

    ~ThreadParams() {
        delete[] mutualInfo;
    }

    const int tid;
    Dataset &dataset;
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    const uint16_t numOutputs;
    MutualInfo *mutualInfo;
    Statistics &statistics;
};


#endif /* THREADPARAMS_H_ */
