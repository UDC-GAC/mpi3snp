//
// Created by chpon on 26/09/2018.
//

#include "CUDAError.h"
#include <cuda_runtime.h>

CUDAError::CUDAError() : Engine::Error(
        std::string(cudaGetErrorName(cudaGetLastError())) + ": " + cudaGetErrorString(cudaGetLastError())) {
}

CUDAError::CUDAError(const std::string &message) : Engine::Error(message) {
}
