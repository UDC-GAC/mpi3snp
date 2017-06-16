//
// Created by christian on 14/06/17.
//

#ifndef MPI3SNP_DATASET_H
#define MPI3SNP_DATASET_H

#include <string>
#include <stdexcept>
#include <vector>
#include "Individual.h"
#include "SNP.h"

class Dataset {
public:
    class ReadError : public std::runtime_error {
    public:
        ReadError(const std::string message) : runtime_error(message) {}

        virtual ~ReadError() {};
    };

    Dataset(std::string tped_path, std::string tfam_path);

    ~Dataset();

    inline std::vector<uint32_t> *&Get_cases(){
        return cases;
    }

    inline std::vector<uint32_t> *&Get_ctrls(){
        return ctrls;
    }

    inline uint32_t Get_SNP_count(){
        return snp_count;
    }

    inline uint16_t Get_ctrl_count(){
        return num_ctrls;
    }

    inline uint16_t Get_case_count(){
        return num_cases;
    }

private:

    void Bitvector_representation(std::vector<Individual> &inds, std::vector<SNP> &snps);

    std::vector<uint32_t> *cases;
    std::vector<uint32_t> *ctrls;
    uint32_t snp_count;
    uint16_t num_cases, num_ctrls;
};

#endif //MPI3SNP_DATASET_H
