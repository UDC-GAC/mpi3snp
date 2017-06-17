/*
 * GPUContTable.h
 *
 *  Created on: 09/10/2014
 *      Author: jorge
 */

#ifndef CONTTABLE_H_
#define CONTTABLE_H_

#include "CUDAError.h"
#include "Macros.h"

/*
 * Structure for the auxiliar 2-SNP contingency tables
 */

struct GPUDoubleContTable {
    GPUDoubleContTable() {

    }

    ~GPUDoubleContTable() {
    }

    void initialize(uint16_t numEntriesCase, uint16_t numEntriesCtrl) {

        if (cudaSuccess != cudaMalloc(&_cases00, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases01, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases02, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases10, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases11, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases12, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases20, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases21, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_cases22, numEntriesCase * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls00, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls01, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls02, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls10, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls11, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls12, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls20, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls21, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());

        if (cudaSuccess != cudaMalloc(&_ctrls22, numEntriesCtrl * sizeof(uint32_t)))
            throw CUDAError(cudaGetLastError());
    }

    void finalize() {
        if (_cases00) {
            if (cudaSuccess != cudaFree(_cases00))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases01) {
            if (cudaSuccess != cudaFree(_cases01))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases02) {
            if (cudaSuccess != cudaFree(_cases02))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases10) {
            if (cudaSuccess != cudaFree(_cases10))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases11) {
            if (cudaSuccess != cudaFree(_cases11))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases12) {
            if (cudaSuccess != cudaFree(_cases12))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases20) {
            if (cudaSuccess != cudaFree(_cases20))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases21) {
            if (cudaSuccess != cudaFree(_cases21))
                throw CUDAError(cudaGetLastError());
        }
        if (_cases22) {
            if (cudaSuccess != cudaFree(_cases22))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls00) {
            if (cudaSuccess != cudaFree(_ctrls00))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls01) {
            if (cudaSuccess != cudaFree(_ctrls01))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls02) {
            if (cudaSuccess != cudaFree(_ctrls02))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls10) {
            if (cudaSuccess != cudaFree(_ctrls10))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls11) {
            if (cudaSuccess != cudaFree(_ctrls11))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls12) {
            if (cudaSuccess != cudaFree(_ctrls12))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls20) {
            if (cudaSuccess != cudaFree(_ctrls20))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls21) {
            if (cudaSuccess != cudaFree(_ctrls21))
                throw CUDAError(cudaGetLastError());
        }
        if (_ctrls22) {
            if (cudaSuccess != cudaFree(_ctrls22))
                throw CUDAError(cudaGetLastError());
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
};


#endif /* CONTTABLE_H_ */
