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

    void run(std::vector<MutualInfo> *mutual_info);

private:
    GPUSNPDistributor *distributor;
    std::vector<uint16_t> gpu_ids;
    uint16_t num_outputs;
    bool is_mi;

    static void *handle(void *arg);
};

#endif /* GPUENGINE_H_ */
