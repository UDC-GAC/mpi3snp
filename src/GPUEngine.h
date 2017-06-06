/*
 * GPUEngine.h
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#include "Options.h"
#include "GPUSNPDistributor.h"

class GPUEngine {
public:
    GPUEngine(Options *options);

    ~GPUEngine();

    void run();

private:
    GPUSNPDistributor *distributor;
    uint16_t _gpuId;
    int num_gpus;
    int *gpu_id;
    int num_outputs;
    bool is_mi;

    static void *handle(void *arg);
};

#endif /* GPUENGINE_H_ */
