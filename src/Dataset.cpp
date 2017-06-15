//
// Created by christian on 14/06/17.
//

#include "Dataset.h"
#include <fstream>
#include <iostream>

Dataset::Dataset(std::string tped_path, std::string tfam_path) :
        cases(new std::vector<uint32_t>[3]),
        ctrls(new std::vector<uint32_t>[3]) {
    std::ifstream file;
    std::vector<Individual> individuals;
    std::vector<AltSNP> snps;

    try {
        file.open(tfam_path.c_str(), std::ios::in);
        file.exceptions(std::ifstream::badbit);
        Individual ind;
        while (file >> ind) {
            individuals.push_back(ind);
        }
    } catch (const std::ifstream::failure &e) {
        throw ReadError(e.what());
    } catch (const Individual::InvalidIndividual &e) {
        throw ReadError("Error in " + tfam_path + ":" + std::to_string(individuals.size() + 1) + ": " + e.what());
    }
    file.close();

    try {
        file.open(tped_path, std::ios::in);
        file.exceptions(std::ifstream::badbit);

        AltSNP snp;
        while (file >> snp) {
            snps.push_back(snp);
        }
    } catch (const std::ifstream::failure &e) {
        throw ReadError(e.what());
    } catch (const AltSNP::InvalidSNP &e) {
        throw ReadError("Error in " + tped_path + ":" + std::to_string(snps.size() + 1) + ": " + e.what());
    }
    file.close();

    Bitvector_representation(individuals, snps);
}

void Dataset::Bitvector_representation(std::vector<Individual> &inds, std::vector<AltSNP> &snps) {
    uint32_t buffer[3];
    int count;
    int i, j, k;

    // Iterate on all SNPs from an specific individual, for all individuals
    for (i = 0; i < inds.size(); i++) {
        j = 0;
        while (j < snps.size()) {
            // Reset buffers
            count = 0;
            for (k = 0; k < 3; k++) {
                buffer[k] = 0;
            }

            // Read all SNPs for a given individual until buffers are full or all SNPs are read
            while (j < snps.size() && count < 32) {
                for (k = 0; k < 3; k++) {
                    buffer[k] = buffer[k] << 1;
                    buffer[k] += snps[j].genotypes[i] == k;
                }
                j++;
                count++;
            }

            // Shift all bits left if 32-bit word is not full
            if (count < 32) {
                for (k = 0; k < 3; k++) {
                    buffer[k] = buffer[k] << (32 - count);
                }
            }

            // Store 32-bit words in each corresponding vector
            if (inds[i].ph == 1) {
                for (k = 0; k < 3; k++) {
                    ctrls[k].push_back(buffer[k]);
                }
            } else {
                for (k = 0; k < 3; k++) {
                    cases[k].push_back(buffer[k]);
                }
            }
        }
    }
}

Dataset::~Dataset() {
    delete[] cases;
    delete[] ctrls;
}

std::vector<uint32_t> *&Dataset::Get_cases() {
    return cases;
}

std::vector<uint32_t> *&Dataset::Get_ctrls() {
    return ctrls;
}