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

        Builder &Set_gpu_ids(std::vector<unsigned int> gpu_ids);

        Builder &Set_cpu_threads(unsigned int threads);

        Builder &Set_use_mi(bool use_mi);

        Search *Create_object();

    private:
        Search *search_obj;
    };

    ~Search();

    // Execute the epistasis search
    void execute();

private:
    Search();

    std::string tped_file, tfam_file, out_file;
    unsigned int num_outputs;
    std::vector<unsigned int> gpu_ids;
    unsigned int cpu_threads;
    bool use_mi;

    MPI_Datatype MPI_MUTUAL_INFO;
};

#endif /* GPUSEARCHMI_H_ */