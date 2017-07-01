//
// Created by christian on 27/06/17.
//

#ifndef MPI3SNP_ENTROPYSEARCH_H
#define MPI3SNP_ENTROPYSEARCH_H

#include <cstdint>
#include <vector>
#include "../MutualInfo.h"
#include "ContTable.h"

class EntropySearch {
public:
    EntropySearch(uint32_t numSNPs, uint16_t numCases, const std::vector<std::vector<uint32_t> *> &cases,
                  uint16_t numCtrls, const std::vector<std::vector<uint32_t> *> &ctrls);

    uint64_t mutualInfo(const std::vector<std::pair<uint32_t, uint32_t>> &pairs, MutualInfo *mutualInfo,
                        uint16_t numOutputs, float &minMI, uint16_t &minMIPos, uint16_t &numEntriesWithMI);

private:

    void _fillDoubleContTable(std::vector<uint32_t> *s1_ctrls, std::vector<uint32_t> *s1_cases,
                              std::vector<uint32_t> *s2_ctrls, std::vector<uint32_t> *s2_cases, DoubleContTable *table);

    float _calcDoubleMI(DoubleContTable *table);

    void
    _fillTripleContTable(DoubleContTable *doubleTable, TripleContTable *tripleTable, std::vector<uint32_t> *s3_ctrls,
                         std::vector<uint32_t> *s3_cases);

    float _calcTripleMI(TripleContTable *table);

    const uint32_t num_snps;
    const std::vector<std::vector<uint32_t> *> &cases;
    const uint16_t num_cases;
    const std::vector<std::vector<uint32_t> *> &ctrls;
    const uint16_t num_ctrls;
    const uint16_t _numEntriesCases;
    const uint16_t _numEntriesCtrls;
    float invInds;
    // Entropy of Y
    float entY;

};


#endif //MPI3SNP_ENTROPYSEARCH_H
