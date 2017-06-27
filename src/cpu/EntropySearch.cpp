//
// Created by christian on 27/06/17.
//

#include "EntropySearch.h"

EntropySearch::EntropySearch(uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls) {
    _numSNPs = numSNPs;
    _numCases = numCases;
    _numCtrls = numCtrls;
    _numEntriesCases = numCases / 32 + ((numCases % 32) > 0);
    _numEntriesCtrls = numCtrls / 32 + ((numCtrls % 32) > 0);
    _invInds = 1.0 / (numCases + numCtrls);

    float p = numCases * _invInds;
    _entY = (-1.0) * p * log2(p);

    p = numCtrls * _invInds;
    _entY -= p * log2(p);
}

uint64_t EntropySearch::mutualInfo(vector<SNP2 *> snpSet, uint64_t numPairs, uint32_t *ids, MutualInfo *mutualInfo,
                                   uint16_t numOutputs, float &minMI, uint16_t &minMIPos, uint16_t &numEntriesWithMI) {

    uint32_t id1, id2, id3;
    float auxMIValue;
    uint64_t numAnal = 0;
    MutualInfo *auxMI;
    DoubleContTable *doubleTable = new DoubleContTable(_numCases, _numCtrls);
    TripleContTable tripleTable;

#ifdef BENCHMARKING_PARTS
    double partTime;
#endif

    for (uint64_t iterPairs = 0; iterPairs < numPairs; iterPairs++) {

        id1 = ids[2 * iterPairs];
        id2 = ids[2 * iterPairs + 1];

        _fillDoubleContTable(snpSet[id1], snpSet[id2], doubleTable);

#ifdef DOUBLE_SNP_ANALYSIS
        auxMIValue = _calcDoubleMI(doubleTable);

        // There are empty values in the array
        if (numEntriesWithMI < numOutputs) {
            auxMI = &mutualInfo[numEntriesWithMI];
            auxMI->_id1 = id1;
            auxMI->_id2 = id2;
            auxMI->_id3 = 0;
            auxMI->_mutualInfoValue = auxMIValue;

            // If this is the minimum value of the array
            if (auxMIValue < minMI) {
                minMI = auxMIValue;
                minMIPos = numEntriesWithMI;
            }

            numEntriesWithMI++;
        } else if (auxMIValue > minMI) { // The value must be inserted
            auxMI = &mutualInfo[minMIPos];
            auxMI->_id1 = id1;
            auxMI->_id2 = id2;
            auxMI->_id3 = 0;
            auxMI->_mutualInfoValue = auxMIValue;

            // Find the new minimum
            auxMI = min_element(mutualInfo, mutualInfo + numOutputs);
            minMI = auxMI->_mutualInfoValue;
            uint16_t i = 0;
            while (1) {
                if (mutualInfo[i]._mutualInfoValue == minMI) {
                    break;
                }
                i++;
            }
            minMIPos = i;
        }
#endif

#ifdef TRIPLE_SNP_ANALYSIS
        for (id3 = id2 + 1; id3 < _numSNPs; id3++) {
            // Generate the contingency table of the 3-SNP
            _fillTripleContTable(doubleTable, &tripleTable, snpSet[id3]);

            // Calculate the MI with the contingency table
            auxMIValue = _calcTripleMI(&tripleTable);

            // There are empty values in the array
            if (numEntriesWithMI < numOutputs) {
                auxMI = &mutualInfo[numEntriesWithMI];
                auxMI->_id1 = id1;
                auxMI->_id2 = id2;
                auxMI->_id3 = id3;
                auxMI->_mutualInfoValue = auxMIValue;

                // If this is the minimum value of the array
                if (auxMIValue < minMI) {
                    minMI = auxMIValue;
                    minMIPos = numEntriesWithMI;
                }

                numEntriesWithMI++;
            } else if (auxMIValue > minMI) { // The value must be inserted
                auxMI = &mutualInfo[minMIPos];
                auxMI->_id1 = id1;
                auxMI->_id2 = id2;
                auxMI->_id3 = id3;
                auxMI->_mutualInfoValue = auxMIValue;

                // Find the new minimum
                auxMI = min_element(mutualInfo, mutualInfo + numOutputs);
                minMI = auxMI->_mutualInfoValue;
                uint16_t i = 0;
                while (1) {
                    if (mutualInfo[i]._mutualInfoValue == minMI) {
                        break;
                    }
                    i++;
                }
                minMIPos = i;
            }

            numAnal++;
        }
#endif
    }

    delete doubleTable;
    return numAnal;
}

void EntropySearch::_fillDoubleContTable(SNP2 *snp1, SNP2 *snp2, DoubleContTable *table) {

    uint16_t sum00 = 0, sum01 = 0, sum02 = 0, sum10 = 0, sum11 = 0, sum12 = 0, sum20 = 0, sum21 = 0, sum22 = 0;

    for (int i = 0; i < _numEntriesCases; i++) {
        table->_cases00[i] = snp1->_case0Values[i] & snp2->_case0Values[i];
        table->_cases01[i] = snp1->_case0Values[i] & snp2->_case1Values[i];
        table->_cases02[i] = snp1->_case0Values[i] & snp2->_case2Values[i];
        table->_cases10[i] = snp1->_case1Values[i] & snp2->_case0Values[i];
        table->_cases11[i] = snp1->_case1Values[i] & snp2->_case1Values[i];
        table->_cases12[i] = snp1->_case1Values[i] & snp2->_case2Values[i];
        table->_cases20[i] = snp1->_case2Values[i] & snp2->_case0Values[i];
        table->_cases21[i] = snp1->_case2Values[i] & snp2->_case1Values[i];
        table->_cases22[i] = snp1->_case2Values[i] & snp2->_case2Values[i];
    }

    for (int i = 0; i < _numEntriesCtrls; i++) {
        table->_ctrls00[i] = snp1->_ctrl0Values[i] & snp2->_ctrl0Values[i];
        table->_ctrls01[i] = snp1->_ctrl0Values[i] & snp2->_ctrl1Values[i];
        table->_ctrls02[i] = snp1->_ctrl0Values[i] & snp2->_ctrl2Values[i];
        table->_ctrls10[i] = snp1->_ctrl1Values[i] & snp2->_ctrl0Values[i];
        table->_ctrls11[i] = snp1->_ctrl1Values[i] & snp2->_ctrl1Values[i];
        table->_ctrls12[i] = snp1->_ctrl1Values[i] & snp2->_ctrl2Values[i];
        table->_ctrls20[i] = snp1->_ctrl2Values[i] & snp2->_ctrl0Values[i];
        table->_ctrls21[i] = snp1->_ctrl2Values[i] & snp2->_ctrl1Values[i];
        table->_ctrls22[i] = snp1->_ctrl2Values[i] & snp2->_ctrl2Values[i];
    }
}

float EntropySearch::_calcDoubleMI(DoubleContTable *table) {

    uint32_t contCases[9];
    uint32_t contCtrls[9];

    for (int i = 0; i < 9; i++) {
        contCases[i] = 0;
        contCtrls[i] = 0;
    }

    for (int i = 0; i < _numEntriesCases; i++) {
        contCases[0] += Utils::popcount(table->_cases00[i]);
        contCases[1] += Utils::popcount(table->_cases01[i]);
        contCases[2] += Utils::popcount(table->_cases02[i]);
        contCases[3] += Utils::popcount(table->_cases10[i]);
        contCases[4] += Utils::popcount(table->_cases11[i]);
        contCases[5] += Utils::popcount(table->_cases12[i]);
        contCases[6] += Utils::popcount(table->_cases20[i]);
        contCases[7] += Utils::popcount(table->_cases21[i]);
        contCases[8] += Utils::popcount(table->_cases22[i]);
    }
    for (int i = 0; i < _numEntriesCtrls; i++) {
        contCtrls[0] += Utils::popcount(table->_ctrls00[i]);
        contCtrls[1] += Utils::popcount(table->_ctrls01[i]);
        contCtrls[2] += Utils::popcount(table->_ctrls02[i]);
        contCtrls[3] += Utils::popcount(table->_ctrls10[i]);
        contCtrls[4] += Utils::popcount(table->_ctrls11[i]);
        contCtrls[5] += Utils::popcount(table->_ctrls12[i]);
        contCtrls[6] += Utils::popcount(table->_ctrls20[i]);
        contCtrls[7] += Utils::popcount(table->_ctrls21[i]);
        contCtrls[8] += Utils::popcount(table->_ctrls22[i]);
    }

    float entX = 0.0;
    float entAll = 0.0;
    float pCase, pCtrl;

    for (int i = 0; i < 9; i++) {
        pCase = contCases[i] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }

        pCtrl = contCtrls[i] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }
    }

    return entX + _entY - entAll;
}

void EntropySearch::_fillTripleContTable(DoubleContTable *doubleTable, TripleContTable *tripleTable, SNP2 *snp3) {

    uint32_t aux;
    uint32_t auxSNP3Value;

    tripleTable->clearValues();

    for (int i = 0; i < _numEntriesCases; i++) {
        auxSNP3Value = snp3->_case0Values[i];

        aux = doubleTable->_cases00[i] & auxSNP3Value;
        tripleTable->_cases[0] += Utils::popcount(aux);

        aux = doubleTable->_cases01[i] & auxSNP3Value;
        tripleTable->_cases[1] += Utils::popcount(aux);

        aux = doubleTable->_cases02[i] & auxSNP3Value;
        tripleTable->_cases[2] += Utils::popcount(aux);

        aux = doubleTable->_cases10[i] & auxSNP3Value;
        tripleTable->_cases[3] += Utils::popcount(aux);

        aux = doubleTable->_cases11[i] & auxSNP3Value;
        tripleTable->_cases[4] += Utils::popcount(aux);

        aux = doubleTable->_cases12[i] & auxSNP3Value;
        tripleTable->_cases[5] += Utils::popcount(aux);

        aux = doubleTable->_cases20[i] & auxSNP3Value;
        tripleTable->_cases[6] += Utils::popcount(aux);

        aux = doubleTable->_cases21[i] & auxSNP3Value;
        tripleTable->_cases[7] += Utils::popcount(aux);

        aux = doubleTable->_cases22[i] & auxSNP3Value;
        tripleTable->_cases[8] += Utils::popcount(aux);


        auxSNP3Value = snp3->_case1Values[i];

        aux = doubleTable->_cases00[i] & auxSNP3Value;
        tripleTable->_cases[9] += Utils::popcount(aux);

        aux = doubleTable->_cases01[i] & auxSNP3Value;
        tripleTable->_cases[10] += Utils::popcount(aux);

        aux = doubleTable->_cases02[i] & auxSNP3Value;
        tripleTable->_cases[11] += Utils::popcount(aux);

        aux = doubleTable->_cases10[i] & auxSNP3Value;
        tripleTable->_cases[12] += Utils::popcount(aux);

        aux = doubleTable->_cases11[i] & auxSNP3Value;
        tripleTable->_cases[13] += Utils::popcount(aux);

        aux = doubleTable->_cases12[i] & auxSNP3Value;
        tripleTable->_cases[14] += Utils::popcount(aux);

        aux = doubleTable->_cases20[i] & auxSNP3Value;
        tripleTable->_cases[15] += Utils::popcount(aux);

        aux = doubleTable->_cases21[i] & auxSNP3Value;
        tripleTable->_cases[16] += Utils::popcount(aux);

        aux = doubleTable->_cases22[i] & auxSNP3Value;
        tripleTable->_cases[17] += Utils::popcount(aux);


        auxSNP3Value = snp3->_case2Values[i];

        aux = doubleTable->_cases00[i] & auxSNP3Value;
        tripleTable->_cases[18] += Utils::popcount(aux);

        aux = doubleTable->_cases01[i] & auxSNP3Value;
        tripleTable->_cases[19] += Utils::popcount(aux);

        aux = doubleTable->_cases02[i] & auxSNP3Value;
        tripleTable->_cases[20] += Utils::popcount(aux);

        aux = doubleTable->_cases10[i] & auxSNP3Value;
        tripleTable->_cases[21] += Utils::popcount(aux);

        aux = doubleTable->_cases11[i] & auxSNP3Value;
        tripleTable->_cases[22] += Utils::popcount(aux);

        aux = doubleTable->_cases12[i] & auxSNP3Value;
        tripleTable->_cases[23] += Utils::popcount(aux);

        aux = doubleTable->_cases20[i] & auxSNP3Value;
        tripleTable->_cases[24] += Utils::popcount(aux);

        aux = doubleTable->_cases21[i] & auxSNP3Value;
        tripleTable->_cases[25] += Utils::popcount(aux);

        aux = doubleTable->_cases22[i] & auxSNP3Value;
        tripleTable->_cases[26] += Utils::popcount(aux);
    }

    for (int i = 0; i < _numEntriesCtrls; i++) {
        auxSNP3Value = snp3->_ctrl0Values[i];

        aux = doubleTable->_ctrls00[i] & auxSNP3Value;
        tripleTable->_ctrls[0] += Utils::popcount(aux);

        aux = doubleTable->_ctrls01[i] & auxSNP3Value;
        tripleTable->_ctrls[1] += Utils::popcount(aux);

        aux = doubleTable->_ctrls02[i] & auxSNP3Value;
        tripleTable->_ctrls[2] += Utils::popcount(aux);

        aux = doubleTable->_ctrls10[i] & auxSNP3Value;
        tripleTable->_ctrls[3] += Utils::popcount(aux);

        aux = doubleTable->_ctrls11[i] & auxSNP3Value;
        tripleTable->_ctrls[4] += Utils::popcount(aux);

        aux = doubleTable->_ctrls12[i] & auxSNP3Value;
        tripleTable->_ctrls[5] += Utils::popcount(aux);

        aux = doubleTable->_ctrls20[i] & auxSNP3Value;
        tripleTable->_ctrls[6] += Utils::popcount(aux);

        aux = doubleTable->_ctrls21[i] & auxSNP3Value;
        tripleTable->_ctrls[7] += Utils::popcount(aux);

        aux = doubleTable->_ctrls22[i] & auxSNP3Value;
        tripleTable->_ctrls[8] += Utils::popcount(aux);


        auxSNP3Value = snp3->_ctrl1Values[i];

        aux = doubleTable->_ctrls00[i] & auxSNP3Value;
        tripleTable->_ctrls[9] += Utils::popcount(aux);

        aux = doubleTable->_ctrls01[i] & auxSNP3Value;
        tripleTable->_ctrls[10] += Utils::popcount(aux);

        aux = doubleTable->_ctrls02[i] & auxSNP3Value;
        tripleTable->_ctrls[11] += Utils::popcount(aux);

        aux = doubleTable->_ctrls10[i] & auxSNP3Value;
        tripleTable->_ctrls[12] += Utils::popcount(aux);

        aux = doubleTable->_ctrls11[i] & auxSNP3Value;
        tripleTable->_ctrls[13] += Utils::popcount(aux);

        aux = doubleTable->_ctrls12[i] & auxSNP3Value;
        tripleTable->_ctrls[14] += Utils::popcount(aux);

        aux = doubleTable->_ctrls20[i] & auxSNP3Value;
        tripleTable->_ctrls[15] += Utils::popcount(aux);

        aux = doubleTable->_ctrls21[i] & auxSNP3Value;
        tripleTable->_ctrls[16] += Utils::popcount(aux);

        aux = doubleTable->_ctrls22[i] & auxSNP3Value;
        tripleTable->_ctrls[17] += Utils::popcount(aux);


        auxSNP3Value = snp3->_ctrl2Values[i];

        aux = doubleTable->_ctrls00[i] & auxSNP3Value;
        tripleTable->_ctrls[18] += Utils::popcount(aux);

        aux = doubleTable->_ctrls01[i] & auxSNP3Value;
        tripleTable->_ctrls[19] += Utils::popcount(aux);

        aux = doubleTable->_ctrls02[i] & auxSNP3Value;
        tripleTable->_ctrls[20] += Utils::popcount(aux);

        aux = doubleTable->_ctrls10[i] & auxSNP3Value;
        tripleTable->_ctrls[21] += Utils::popcount(aux);

        aux = doubleTable->_ctrls11[i] & auxSNP3Value;
        tripleTable->_ctrls[22] += Utils::popcount(aux);

        aux = doubleTable->_ctrls12[i] & auxSNP3Value;
        tripleTable->_ctrls[23] += Utils::popcount(aux);

        aux = doubleTable->_ctrls20[i] & auxSNP3Value;
        tripleTable->_ctrls[24] += Utils::popcount(aux);

        aux = doubleTable->_ctrls21[i] & auxSNP3Value;
        tripleTable->_ctrls[25] += Utils::popcount(aux);

        aux = doubleTable->_ctrls22[i] & auxSNP3Value;
        tripleTable->_ctrls[26] += Utils::popcount(aux);
    }
}

float EntropySearch::_calcTripleMI(TripleContTable *table) {

    float entX = 0.0;
    float entAll = 0.0;
    float pCase, pCtrl;

    for (int i = 0; i < 27; i++) {

        pCase = table->_cases[i] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }

        pCtrl = table->_ctrls[i] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }
    }

    return entX + _entY - entAll;
}