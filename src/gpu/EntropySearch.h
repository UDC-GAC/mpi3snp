//
// Created by christian on 06/06/17.
//

#ifndef MPI3SNP_ENTROPYSEARCH_H
#define MPI3SNP_ENTROPYSEARCH_H

#include "../MutualInfo.h"
#include "GPUContTable.h"
#include <vector>

class EntropySearch {
public:
    EntropySearch(bool isMI, uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls, uint16_t numOutputs,
                  unsigned long block_size, std::vector<std::vector<uint32_t> *> cases,
                  std::vector<std::vector<uint32_t> *> ctrls);

    ~EntropySearch();

    void mutualInfo(uint64_t numPairs, uint2 *ids, MutualInfo *mutualInfo, float &minMI, uint16_t &minMIPos,
                    uint16_t &numEntriesWithMI);

private:
    bool _isMI;
    const unsigned long block_size;
    uint16_t _numOutputs;

    uint32_t _numSNPs;
    uint16_t _numCases;
    uint16_t _numCtrls;
    uint16_t _numEntriesCase;
    uint16_t _numEntriesCtrl;

    uint32_t *_dev0Cases;
    uint32_t *_dev1Cases;
    uint32_t *_dev2Cases;
    uint32_t *_dev0Ctrls;
    uint32_t *_dev1Ctrls;
    uint32_t *_dev2Ctrls;
    uint2 *_devIds;

    // Auxiliary array for the contingency tables between the two kernels
    GPUDoubleContTable *_devDoubleTables;
    GPUDoubleContTable *_tables;

    // Auxiliary array to store the MI values of each block
    float *_devMIValues;
    float *_hostMIValues;

    // Auxiliary arrays to store the ids that are in the list of MIs
    uint3 *_devMiIds;
    uint3 *_hostMiIds;

    void _findNHighestMI(MutualInfo *mutualInfo, uint64_t totalValues, float &minMI, uint16_t &minMIPos,
                         uint16_t &numEntriesWithMI);
};


#endif //MPI3SNP_ENTROPYSEARCH_H
