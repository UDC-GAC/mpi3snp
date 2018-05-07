/*
 * This file is part of MPI3SNP.
 * Copyright (C) 2018 by Christian Ponte
 *
 * MPI3SNP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MPI3SNP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file Dataset.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Dataset class declaration, responsible for reading the data and storing it using a bitwise representation.
 */

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

    enum Representation {
        Regular,
        Transposed
    };

    Dataset(std::string tped_path, std::string tfam_path, Representation rep);

    ~Dataset();

    inline std::vector<std::vector<uint32_t> *> &Get_cases() {
        return cases;
    }

    inline std::vector<std::vector<uint32_t> *> &Get_ctrls() {
        return ctrls;
    }

    inline uint32_t Get_SNP_count() {
        return snp_count;
    }

    inline uint16_t Get_ctrl_count() {
        return num_ctrls;
    }

    inline uint16_t Get_case_count() {
        return num_cases;
    }

private:
    void Regular_representation(std::vector<Individual> &inds, std::vector<SNP> &snps);

    void Transposed_representation(std::vector<Individual> &inds, std::vector<SNP> &snps);

    std::vector<std::vector<uint32_t> *> cases;
    std::vector<std::vector<uint32_t> *> ctrls;
    uint32_t snp_count;
    uint16_t num_cases, num_ctrls;
};

#endif //MPI3SNP_DATASET_H
