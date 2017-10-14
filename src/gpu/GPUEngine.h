/*
 * GPUEngine.h
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#include "../Engine.h"

class GPUEngine : public Engine {
public:
    GPUEngine(unsigned int proc_num, unsigned int proc_id, std::vector<std::pair<unsigned int, unsigned int>> gpu_map,
              bool use_mi);

    void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs,
             Statistics &statistics) override;

private:
    unsigned int proc_num, proc_id, gpu_id;
    bool use_mi;
};

#endif /* GPUENGINE_H_ */
