//
// Created by christian on 17/06/17.
//

#ifndef MPI3SNP_CUDAERROR_H
#define MPI3SNP_CUDAERROR_H

#include <cuda_runtime.h>
#include "../Engine.h"

class CUDAError : public Engine::Error {
public:
    CUDAError(const std::string &message) :
            Engine::Error(message) {};

    CUDAError(const cudaError_t &code) :
            Engine::Error(std::string(cudaGetErrorName(code)) + ": " + cudaGetErrorString(code)) {};

    ~CUDAError() override = default;
};

#endif //MPI3SNP_CUDAERROR_H
