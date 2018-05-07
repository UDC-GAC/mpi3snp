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
 * @file Arg_parser.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Argument parser class. Takes the command line arguments as the input, and returns a structure containing the
 * configuration to be used.
 */

#ifndef MPI3SNP_ARG_PARSER_H
#define MPI3SNP_ARG_PARSER_H

#include <map>
#include <vector>
#include <thread>

class Arg_parser {
public:
    struct Arguments {
        std::string tped;
        std::string tfam;
        std::string output;
        std::vector<std::pair<unsigned int, unsigned int>> gpu_map;
        unsigned int cpu_threads = std::thread::hardware_concurrency();
        bool use_mi = true;
        unsigned int output_num = 10;
        bool debug;
        bool benchmarking;
    };

    class Finalize : public std::runtime_error {
    public:
        explicit Finalize(const std::string &message) : runtime_error(message) {};

        ~Finalize() override = default;
    };

    Arg_parser(int &argc, char **argv);

    Arguments get_arguments();

private:
    const int argc;
    char **argv;
};

#endif //MPI3SNP_ARG_PARSER_H
