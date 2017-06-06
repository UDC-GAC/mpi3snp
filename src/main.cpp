#include "IOMpi.h"
#include "GPUSearchMI.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    double stime, etime;

    /*get the startup time*/
    stime = MPI_Wtime();

    Options options;

    /*parse the arguments*/
    if (!options.parse(argc, argv)) {
        options.printUsage();
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    Search *search = new GPUSearchMI(&options);

    search->execute();

    etime = MPI_Wtime();
    IOMpi::Instance().Mprintf("Overall time: %.2f seconds\n", etime - stime);

    delete search;

    IOMpi::Instance().Deallocate_MPI_resources();
    MPI_Finalize();
    return 0;
}
