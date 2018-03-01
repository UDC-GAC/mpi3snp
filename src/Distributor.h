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
 * @brief Distributor class declaration. This class implements the parallel distribution strategy, shared by both cpu
 * and gpu implementation.
 */

#ifndef MPI3SNP_DISTRIBUTOR_H
#define MPI3SNP_DISTRIBUTOR_H

#include <inttypes.h>
#include <vector>

class Distributor {
public:
    Distributor(unsigned int p_num, unsigned int p_id, unsigned long count);

    void Get_pairs(unsigned int t_num, unsigned int t_id,
                            std::vector<std::pair<uint32_t, uint32_t>> &values);

private:
    const unsigned long count;
    std::vector<uint32_t> distribution;
};

#endif //MPI3SNP_DISTRIBUTOR_H
