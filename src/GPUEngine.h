/*
 * GPUEngine.h
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#include "GPUSNPDistributor.h"
#include "Statistics.h"

class GPUEngine {
public:
    GPUEngine(std::string tped, std::string tfam, int proc_num, int proc_id, std::vector<unsigned int> gpu_ids,
              uint16_t num_outputs, bool use_mi);

    ~GPUEngine();

    void run(std::vector<MutualInfo> &mutual_info, Statistics &statistics);

private:
    GPUSNPDistributor *distributor;
    std::vector<unsigned int> gpu_ids;
    uint16_t num_outputs;
    bool use_mi;

    static void *handle(void *arg);
};

#endif /* GPUENGINE_H_ */
