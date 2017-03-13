#include "SearchMI.h"
#include "IOMpi.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    double stime, etime;

    /*get the startup time*/
    stime = MPI_Wtime();

    Options *options = new Options(&argc, &argv);

    /*parse the arguments*/
    if (!options->parse()) {
        options->printUsage();
        return 0;
    }

    Search *search = new SearchMI(options);

    search->execute();

    MPI_Barrier(MPI_COMM_WORLD);
    etime = MPI_Wtime();

    IOMpi::Instance().Mprintf("Overall time: %.2f seconds\n", etime - stime);

    delete search;
    delete options;

    IOMpi::Deallocate_MPI_resources();
    MPI_Finalize();

    return 0;
}
