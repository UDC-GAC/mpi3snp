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
    EntropySearch(bool isMI, uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls,
                  std::vector<std::vector<uint32_t> *> cases, std::vector<std::vector<uint32_t> *> ctrls);

    ~EntropySearch();

    void mutualInfo(std::vector<std::pair<uint32_t, uint32_t >> pairs, size_t num_outputs, MutualInfo *mutualInfo);

private:
    bool _isMI;

    uint32_t _numSNPs;
    uint16_t _numEntriesCase;
    uint16_t _numEntriesCtrl;

    uint32_t *_dev0Cases;
    uint32_t *_dev1Cases;
    uint32_t *_dev2Cases;
    uint32_t *_dev0Ctrls;
    uint32_t *_dev1Ctrls;
    uint32_t *_dev2Ctrls;
    uint2 *_devIds;

    void _findNHighestMI(uint3 *_hostMiIds, float *_hostMIValues, uint64_t totalValues, float &minMI, uint16_t &minMIPos, uint16_t &numEntriesWithMI,
                         size_t num_outputs, MutualInfo *mutualInfo);
};


#endif //MPI3SNP_ENTROPYSEARCH_H
