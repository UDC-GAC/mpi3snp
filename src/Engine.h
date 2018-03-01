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
 * @file Engine.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Abstract class Engine definition.
 */

#ifndef MPI3SNP_ENGINE_H
#define MPI3SNP_ENGINE_H

#include <vector>
#include "MutualInfo.h"
#include "Statistics.h"

class Engine {
public:
    class Error : public std::runtime_error {
    public:
        explicit Error(const std::string &message) : runtime_error(message) {};

        ~Error() override = default;
    };

    virtual void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs,
                     Statistics &statistics) = 0;
};

#endif //MPI3SNP_ENGINE_H
