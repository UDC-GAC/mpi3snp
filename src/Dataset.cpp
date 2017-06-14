//
// Created by christian on 14/06/17.
//

#include "Dataset.h"
#include <fstream>

Dataset::Dataset(std::string tped_path, std::string tfam_path) {
    std::ifstream file;

    try {
        file.open(tfam_path, std::ios::in);
        file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        Individual ind;
        while (file >> ind) {
            individuals.push_back(ind);
        }
    } catch (const std::ifstream::failure &e) {
        throw ReadError(e.what());
    }
    file.close();

    try {
        file.open(tped_path, std::ios::in);
        file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        SNP snp;
        while (file >> snp) {
            snps.push_back(snp);
        }
    } catch (const std::ifstream::failure &e) {
        throw ReadError(e.what());
    }
    file.close();
}

std::vector<SNP> &Dataset::Get_SNPs() {
    return snps;
}

std::vector<Individual> &Dataset::Get_individuals() {
    return individuals;
}