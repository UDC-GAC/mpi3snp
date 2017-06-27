//
// Created by christian on 27/06/17.
//

#ifndef MPI3SNP_ENTROPYSEARCH_H
#define MPI3SNP_ENTROPYSEARCH_H

#include <cstdint>
#include "SNP2.h"
#include "../MutualInfo.h"
#include "ContTable.h"

class EntropySearch {
public:
    EntropySearch(uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls);

    uint64_t mutualInfo(vector<SNP2 *> snpSet, uint64_t numPairs, uint32_t *ids, MutualInfo *mutualInfo,
                        uint16_t numOutputs, float &minMI, uint16_t &minMIPos, uint16_t &numEntriesWithMI);

private:

    void _fillDoubleContTable(SNP2 *snp1, SNP2 *snp2, DoubleContTable *table);

    float _calcDoubleMI(DoubleContTable *table);

    void _fillTripleContTable(DoubleContTable *doubleTable, TripleContTable *tripleTable, SNP2 *snp3);

    float _calcTripleMI(TripleContTable *table);

    uint32_t _numSNPs;
    uint16_t _numCases;
    uint16_t _numCtrls;
    uint16_t _numEntriesCases;
    uint16_t _numEntriesCtrls;
    float _invInds;
    // Entropy of Y
    float _entY;

};


#endif //MPI3SNP_ENTROPYSEARCH_H
