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

class GPUSearchMI {
public:
    class Builder {
    public:
        Builder(std::string tped_file, std::string tfam_file, std::string out_file);

        Builder &Set_num_outputs(unsigned int num_outputs);

        Builder &Set_gpu_ids(std::vector<unsigned int> gpu_ids);

        Builder &Set_use_mi(bool use_mi);

        GPUSearchMI *Create_object();

    private:
        GPUSearchMI *search_obj;
    };

    ~GPUSearchMI();

    // Execute the epistasis search
    void execute();

private:
    GPUSearchMI();

    std::string tped_file, tfam_file, out_file;
    unsigned int num_outputs;
    std::vector<unsigned int> gpu_ids;
    bool use_mi;

    MPI_Datatype MPI_MUTUAL_INFO;
};

#endif /* GPUSEARCHMI_H_ */
