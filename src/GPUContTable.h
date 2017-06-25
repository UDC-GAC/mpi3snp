/*
 * GPUContTable.h
 *
 *  Created on: 09/10/2014
 *      Author: jorge
 */

#ifndef CONTTABLE_H_
#define CONTTABLE_H_

#include "CUDAError.h"

/*
 * Structure for the auxiliar 2-SNP contingency tables
 */

struct GPUDoubleContTable {
    void initialize(uint16_t numEntriesCase, uint16_t numEntriesCtrl) {
        if (cudaSuccess != cudaMalloc(&base_address, (9 * numEntriesCase + 9 * numEntriesCtrl) * sizeof(uint32_t))) {
            throw CUDAError(cudaGetLastError());
        }

        _cases00 = base_address;
        _cases01 = base_address + numEntriesCase;
        _cases02 = base_address + 2 * numEntriesCase;
        _cases10 = base_address + 3 * numEntriesCase;
        _cases11 = base_address + 4 * numEntriesCase;
        _cases12 = base_address + 5 * numEntriesCase;
        _cases20 = base_address + 6 * numEntriesCase;
        _cases21 = base_address + 7 * numEntriesCase;
        _cases22 = base_address + 8 * numEntriesCase;
        uint32_t *temp_addr = base_address + 9 * numEntriesCase;
        _ctrls00 = temp_addr;
        _ctrls01 = temp_addr + numEntriesCtrl;
        _ctrls02 = temp_addr + 2 * numEntriesCtrl;
        _ctrls10 = temp_addr + 3 * numEntriesCtrl;
        _ctrls11 = temp_addr + 4 * numEntriesCtrl;
        _ctrls12 = temp_addr + 5 * numEntriesCtrl;
        _ctrls20 = temp_addr + 6 * numEntriesCtrl;
        _ctrls21 = temp_addr + 7 * numEntriesCtrl;
        _ctrls22 = temp_addr + 8 * numEntriesCtrl;
    }

    void finalize() {
        if (base_address != nullptr) {
            if (cudaSuccess != cudaFree(base_address)) {
                throw CUDAError(cudaGetLastError());
            }
        }
    }

    uint32_t *_cases00;
    uint32_t *_cases01;
    uint32_t *_cases02;
    uint32_t *_cases10;
    uint32_t *_cases11;
    uint32_t *_cases12;
    uint32_t *_cases20;
    uint32_t *_cases21;
    uint32_t *_cases22;
    uint32_t *_ctrls00;
    uint32_t *_ctrls01;
    uint32_t *_ctrls02;
    uint32_t *_ctrls10;
    uint32_t *_ctrls11;
    uint32_t *_ctrls12;
    uint32_t *_ctrls20;
    uint32_t *_ctrls21;
    uint32_t *_ctrls22;

private:
    uint32_t *base_address = nullptr;
};


#endif /* CONTTABLE_H_ */
