#include "IOMpi.h"
#include "Search.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double stime, etime;

    /*get the startup time*/
    stime = MPI_Wtime();

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
    }

    Search *search = Search::Builder::build_from_args(arguments);
    if (search == nullptr) {
        MPI_Finalize();
        return 0;
    }
    search->execute();

    etime = MPI_Wtime();
    IOMpi::Instance().mprint<IOMpi::B>("Overall time: " + std::to_string(etime - stime) + " seconds\n");

    delete search;

    MPI_Finalize();
    return 0;
}
