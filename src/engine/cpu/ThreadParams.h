/*
 * This file is part of MPI3SNP.
 * Copyright (C) 2018 by Christian Ponte
 *
 * MPI3SNP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MPI3SNP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MPI3SNP_THREADPARAMS_H
#define MPI3SNP_THREADPARAMS_H

#include "MutualInfo.h"
#include "Distributor.h"
#include "Dataset.h"

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


#endif //MPI3SNP_THREADPARAMS_H
