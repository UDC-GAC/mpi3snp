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
 * @brief Distributor class members implementation.
 */

#include "Distributor.h"
#include <cstdlib>

Distributor::Distributor(unsigned int p_num, unsigned int p_id, unsigned long count) :
        count(count) {
    unsigned long dist_size = count / p_num + (p_id < (count % p_num));
    distribution.resize(dist_size);
    distribution[0] = p_id;
    distribution[1] = 2 * p_num - p_id - 1;
    for (unsigned int i = 2; i < dist_size; i++) {
        distribution[i] = distribution[i - 2] + 2 * p_num;
    }
}

void Distributor::Get_pairs(unsigned int t_num, unsigned int t_id,
                                     std::vector<std::pair<uint32_t, uint32_t>> &values) {
    for (size_t i = t_id; i < distribution.size(); i += t_num) {
        for (uint32_t j = distribution[i] + 1; j < count - 1; j++) {
            values.push_back(std::make_pair(distribution[i], j));
        }
    }
}