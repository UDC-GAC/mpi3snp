/*
 * Engine.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include "Macros.h"
#include "SNP2.h"
#include "../MutualInfo.h"
#include "ContTable.h"
#include "SNPDistributor.h"
#include "../Statistics.h"

class CPUEngine {
public:
    CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi);

    virtual ~CPUEngine();

    void execute(std::string tped_file, std::string tfam_file, std::vector<MutualInfo> &mutual_info,
                 uint16_t num_outputs, Statistics &statistics);

private:

    static void *threadMI(void * arg);

    int num_proc;
    int proc_id;
    int num_threads;
    bool use_mi;
};

#endif /* ENGINE_H_ */
