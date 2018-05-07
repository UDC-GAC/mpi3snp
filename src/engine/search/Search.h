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
 * @file Search.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Search class definition.
 */

#ifndef MPI3SNP_SEARCH_H
#define MPI3SNP_SEARCH_H

#include <vector>
#include "Arg_parser.h"
#include "Engine.h"

class Search {
public:
    class Builder {
    public:
        static Search *build_from_args(Arg_parser::Arguments arguments, Statistics &statistics);

        Builder() = delete;
    };

    // Execute the epistasis search
    void execute();

private:
    Search() = default;

    Engine *engine;
    int proc_id, num_proc;
    std::string tped_file, tfam_file, out_file;
    unsigned int num_outputs;
};

#endif //MPI3SNP_SEARCH_H
