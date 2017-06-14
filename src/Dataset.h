//
// Created by christian on 14/06/17.
//

#ifndef MPI3SNP_DATASET_H
#define MPI3SNP_DATASET_H

#include <string>
#include <stdexcept>

class Dataset {
public:
    class ReadError : public std::runtime_error {
    public:
        ReadError(const std::string message) : runtime_error(message) {}

        virtual ~ReadError() {};
    };

    Dataset(std::string tped_path, std::string tfam_path);

    std::vector<SNP> &Get_SNPs();

    std::vector<Individual> &Get_individuals();

private:
    std::vector<SNP> snps;
    std::vector<Individual> individuals;
};

#endif //MPI3SNP_DATASET_H
