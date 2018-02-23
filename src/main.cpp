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
            if (i == info_list.size() - 1){
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
    Search *search = Search::Builder::build_from_args(arguments);
    if (search == nullptr) {
        MPI_Finalize();
        return 0;
    }
    search->execute();

    /*get the ending time*/
    etime = MPI_Wtime();

    IOMpi::Instance().mprint<IOMpi::B>("Overall time: " + std::to_string(etime - stime) + " seconds\n");

    delete search;

    MPI_Finalize();
    return 0;
}
