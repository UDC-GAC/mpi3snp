//
// Created by christian on 10/03/17.
//

#ifndef MPI3SNP_IOMpi_H
#define MPI3SNP_IOMpi_H

#include <mpi.h>

/*
 * CLASS: IOMpi
 *
 * DESCRIPTION: Handles all the I/O throughout the program. Uses Meyers' Singleton pattern, which describes a thread-safe
 *  singleton initialization (C++11 or superior). The object destruction is handled at the end of the program; however,
 *  as MPI dependant resources need to be freed before calling MPI_Finalize, some attributes are destroyed with the call
 *  Deallocate_MPI_resources().
 */

class IOMpi {
public:
    static IOMpi &Instance() {
        static IOMpi inst;
        return inst;
    }

    static void Deallocate_MPI_resources(){
        MPI_Comm_free(&IOMpi::Instance().io_comm);
    }

    // delete copy and move constructors and assign operators
    IOMpi(IOMpi const &) = delete;             // Copy construct
    IOMpi(IOMpi &&) = delete;                  // Move construct
    IOMpi &operator=(IOMpi const &) = delete;  // Copy assign
    IOMpi &operator=(IOMpi &&) = delete;      // Move assign

    int Cprintf(char *format, ...);

protected:
    /* Methods */
    IOMpi();

    ~IOMpi();

    int Get_io_rank();

    /* Attributes */
    const int BUFFER_SIZE = 10240;

    MPI_Comm io_comm;
    int io_rank, my_rank, comm_size;
    int cprintf_tag;
    pthread_mutex_t cprintf_mutex;
    char *io_buff;
};


#endif //MPI3SNP_IOMpi_H
