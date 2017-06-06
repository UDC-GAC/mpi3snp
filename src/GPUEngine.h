/*
 * GPUEngine.h
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#include "Options.h"
#include "GPUContTable.h"
#include "MutualInfo.h"
#include "float.h"
#include "EntropySearch.h"

class GPUEngine {
public:
    inline GPUEngine(Options *options, int gpuId, uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls,
                     uint32_t *host0Cases, uint32_t *host1Cases, uint32_t *host2Cases,
                     uint32_t *host0Ctrls, uint32_t *host1Ctrls, uint32_t *host2Ctrls) {
        _gpuInfo = GPUInfo::getGPUInfo();
        _gpuId = options->getGPUId(gpuId);
        // Set GPU device
        _gpuInfo->setDevice(_gpuId);

        entropySearch = new EntropySearch(options->isMI(), numSNPs, numCases, numCtrls, options->getNumOutputs(),
                                          host0Cases, host1Cases, host2Cases, host0Ctrls, host1Ctrls, host2Ctrls);
    }

    inline ~GPUEngine() {
        delete entropySearch;
    }

    inline void mutualInfo(uint64_t numPairs, uint2 *ids, MutualInfo *mutualInfo, float &minMI, uint16_t &minMIPos,
                           uint16_t &numEntriesWithMI) {
        entropySearch->mutualInfo(numPairs, ids, mutualInfo, minMI, minMIPos, numEntriesWithMI);
    }

private:
    GPUInfo *_gpuInfo;
    uint16_t _gpuId;
    EntropySearch *entropySearch;
};

#endif /* GPUENGINE_H_ */
