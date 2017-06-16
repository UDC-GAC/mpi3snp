//
// Created by christian on 14/06/17.
//

#include "Dataset.h"
#include <fstream>
#include <functional>

Dataset::Dataset(std::string tped_path, std::string tfam_path) {
    std::ifstream file;
    std::vector<Individual> individuals;
    std::vector<SNP> snps;

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

        SNP snp;
        while (file >> snp) {
            if (snp.genotypes.size() == individuals.size()) {
                snps.push_back(snp);
            } else {
                throw ReadError("Error in " + tped_path + ":" + std::to_string(snps.size() + 1) +
                                ": the number of nucleotides does not match the number of individuals");
            }
        }
    } catch (const std::ifstream::failure &e) {
        throw ReadError(e.what());
    } catch (const SNP::InvalidSNP &e) {
        throw ReadError("Error in " + tped_path + ":" + std::to_string(snps.size() + 1) + ": " + e.what());
    }
    file.close();

    cases = new std::vector<uint32_t>[3];
    ctrls = new std::vector<uint32_t>[3];
    snp_count = snps.size();
    Bitvector_representation(individuals, snps);
}

unsigned long find_index(unsigned long start, unsigned long end, std::function<bool(unsigned long)> fun) {
    while (start != end && !fun(start)) {
        start++;
    }
    return start;
}

void Dataset::Bitvector_representation(std::vector<Individual> &inds, std::vector<SNP> &snps) {
    std::vector<unsigned long> scases(32), sctrls(32);
    unsigned long ctrlpos = 0, casepos = 0;
    uint32_t cases_buffer[3], ctrls_buffer[3];
    int i, j;

    num_cases = 0;
    num_ctrls = 0;

    // Iterate on all SNPs considering a block of 32 cases and controls, for all individuals
    while (num_cases + num_ctrls < inds.size()) {
        // Select next 32 controls and cases
        scases.clear();
        while (scases.size() < 32 &&
               (casepos = find_index(casepos, inds.size(), [&inds](int k) { return inds[k].ph == 2; })) < inds.size()) {
            scases.push_back(casepos++);
        }
        sctrls.clear();
        while (sctrls.size() < 32 &&
               (ctrlpos = find_index(ctrlpos, inds.size(), [&inds](int k) { return inds[k].ph == 1; })) < inds.size()) {
            sctrls.push_back(ctrlpos++);
        }
        // Read all SNPs for those 32 controls and cases
        for (i = 0; i < snps.size(); i++) {
            for (j = 0; j < 3; j++) {
                cases_buffer[j] = 0;
                ctrls_buffer[j] = 0;
            }
            for (unsigned long pos : scases) {
                for (j = 0; j < 3; j++) {
                    cases_buffer[j] = cases_buffer[j] << 1;
                    cases_buffer[j] += snps[i].genotypes[pos] == j;
                }
            }
            for (unsigned long pos : sctrls) {
                for (j = 0; j < 3; j++) {
                    ctrls_buffer[j] = ctrls_buffer[j] << 1;
                    ctrls_buffer[j] += snps[i].genotypes[pos] == j;
                }
            }
            // Save buffers when not empty
            if (!scases.empty()) {
                for (j = 0; j < 3; j++) {
                    cases[j].push_back(cases_buffer[j]);
                }
            }
            if (!sctrls.empty()) {
                for (j = 0; j < 3; j++) {
                    ctrls[j].push_back(ctrls_buffer[j]);
                }
            }
        }
        num_cases += scases.size();
        num_ctrls += sctrls.size();
    }
}

Dataset::~Dataset() {
    delete[] cases;
    delete[] ctrls;
}