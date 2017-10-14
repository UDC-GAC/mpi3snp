/*
 * GPUSearchMI.h
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#ifndef GPUSEARCHMI_H_
#define GPUSEARCHMI_H_

#include <mpi.h>
#include <vector>
#include "Arg_parser.h"
#include "Engine.h"

class Search {
public:
    class Builder {
    public:
        static Search *build_from_args(Arg_parser::Arguments arguments);

        Builder() = delete;
    };

    // Execute the epistasis search
    void execute();

private:
    Search() = default;

    Engine *engine;
    int proc_id, num_proc;
    std::string tped_file, tfam_file, out_file;
    unsigned int num_outputs;
};

#endif /* GPUSEARCHMI_H_ */
