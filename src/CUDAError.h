//
// Created by christian on 17/06/17.
//

#ifndef MPI3SNP_CUDAERROR_H
#define MPI3SNP_CUDAERROR_H

#include <cuda_runtime.h>
#include <stdexcept>

class CUDAError : public std::runtime_error {
public:
    CUDAError(const std::string &message) :
            std::runtime_error(message) {};

    CUDAError(const cudaError_t &code) :
            std::runtime_error(std::string(cudaGetErrorName(code)) + ": " + cudaGetErrorString(code)) {};
};

#endif //MPI3SNP_CUDAERROR_H
