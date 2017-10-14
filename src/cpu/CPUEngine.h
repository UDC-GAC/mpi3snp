/*
 * Engine.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include "../Engine.h"

class CPUEngine : public Engine {
public:
    CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi);

    virtual ~CPUEngine() = default;

    void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs,
                  Statistics &statistics) override;

private:
    static void *threadMI(void *arg);

    int num_proc;
    int proc_id;
    int num_threads;
    bool use_mi;
};

#endif /* ENGINE_H_ */
