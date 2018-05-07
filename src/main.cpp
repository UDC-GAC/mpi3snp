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
 * @file main.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Main program. Initializes the MPI environment, parses the command-line arguments, prints debug information if
 * necessary, instantiates the Search and terminates the execution.
 */

#include <mpi.h>
#include "IOMpi.h"
#include "Search.h"
#include "Node_information.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double stime, etime;

    Arg_parser::Arguments arguments;
    try {
        arguments = Arg_parser(argc, argv).get_arguments();
    } catch (const Arg_parser::Finalize &finalize) {
        MPI_Finalize();
        return 0;
    }

    // Set adequate information printing level
    if (arguments.benchmarking) {
        IOMpi::Instance().Set_print_level(IOMpi::B);
    } else if (arguments.debug) {
        IOMpi::Instance().Set_print_level(IOMpi::D);

        // Print node information
        auto info_list = Node_information::gather(0);
        IOMpi::Instance().mprint<IOMpi::D>("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
        for (int i = 0; i < info_list.size(); i++) {
            IOMpi::Instance().mprint<IOMpi::D>(info_list[i]->to_string());
            if (i == info_list.size() - 1) {
                IOMpi::Instance().mprint<IOMpi::D>("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n\n");
            } else {
                IOMpi::Instance().mprint<IOMpi::D>("---------------------------------------------\n");
            }
            delete info_list[i];
        }
        info_list.clear();
    }

    /*get the startup time*/
    stime = MPI_Wtime();

    // Execute search
    Statistics statistics;
    Search *search = Search::Builder::build_from_args(arguments, statistics);
    if (search == nullptr) {
        MPI_Finalize();
        return 0;
    }
    search->execute();

    /*get the ending time*/
    etime = MPI_Wtime();

    delete search;
    // Print runtime statistics to stdout
    IOMpi::Instance().cprint<IOMpi::B>(statistics.To_string());
    IOMpi::Instance().mprint<IOMpi::B>("Overall time: " + std::to_string(etime - stime) + " seconds\n");

    MPI_Finalize();
    return 0;
}
