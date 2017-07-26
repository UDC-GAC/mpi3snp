/*
 * GPUEngine.h
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#include "../Statistics.h"
#include "../MutualInfo.h"

class GPUEngine {
public:
    GPUEngine(unsigned int proc_num, unsigned int proc_id, bool use_mi);

    void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, unsigned int num_outputs,
             Statistics &statistics);

private:
    unsigned int proc_num, proc_id, gpu_id;
    bool use_mi;
};

#endif /* GPUENGINE_H_ */
