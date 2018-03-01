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
 * @file MutualInfo.h
 * @author Christian Ponte
 * @author Jorge González
 * @date 1 March 2018
 *
 * @brief MutualInfo structure definition, used for the mutual information outputs with the result and the SNPs
 * associated.
 */

#ifndef MPI3SNP_MUTUALINFO_H
#define MPI3SNP_MUTUALINFO_H

#include <inttypes.h>
#include <string>

struct MutualInfo {
    MutualInfo() {
        _id1 = 0;
        _id2 = 0;
        _id3 = 0;
        _mutualInfoValue = 0.0;
    }

    uint32_t _id1;
    uint32_t _id2;
    uint32_t _id3;
    float _mutualInfoValue;

    bool operator<(const MutualInfo mi) const { return _mutualInfoValue < mi._mutualInfoValue; }

    inline std::string To_string() {
        return std::to_string(_id1) + " " + std::to_string(_id2) + " " + std::to_string(_id3) + " " +
               std::to_string(_mutualInfoValue);
    }
};

#endif //MPI3SNP_MUTUALINFO_H
