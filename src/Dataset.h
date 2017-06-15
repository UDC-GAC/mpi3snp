//
// Created by christian on 14/06/17.
//

#ifndef MPI3SNP_DATASET_H
#define MPI3SNP_DATASET_H

#include <string>
#include <stdexcept>
#include <vector>
#include "Individual.h"
#include "AltSNP.h"

class Dataset {
public:
    class ReadError : public std::runtime_error {
    public:
        ReadError(const std::string message) : runtime_error(message) {}

        virtual ~ReadError() {};
    };

    Dataset(std::string tped_path, std::string tfam_path);

    ~Dataset();

    std::vector<uint32_t> *&Get_cases();

    std::vector<uint32_t> *&Get_ctrls();

private:

    void Bitvector_representation(std::vector<Individual> &inds, std::vector<AltSNP> &snps);

    std::vector<uint32_t> *cases;
    std::vector<uint32_t> *ctrls;
};

#endif //MPI3SNP_DATASET_H
