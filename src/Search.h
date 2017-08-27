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

class Search {
public:
    class Builder {
    public:
        Builder(std::string tped_file, std::string tfam_file, std::string out_file);

        Builder &Set_num_outputs(unsigned int num_outputs);

        Builder &Set_gpu_map(std::vector<std::pair<unsigned int, unsigned int>> gpu_map);

        Builder &Set_cpu_threads(unsigned int threads);

        Builder &Set_use_mi(bool use_mi);

        Search *Create_object();

    private:
        Search *search_obj;
    };

    // Execute the epistasis search
    void execute();

private:
    Search();

    std::string tped_file, tfam_file, out_file;
    unsigned int num_outputs;
    std::vector<std::pair<unsigned int, unsigned int>> gpu_map;
    unsigned int cpu_threads;
    bool use_mi;
};

#endif /* GPUSEARCHMI_H_ */
