//
// Created by christian on 10/03/17.
//

#ifndef MPI3SNP_IOMpi_H
#define MPI3SNP_IOMpi_H

#include <mpi.h>

class IOMpi {
public:
    IOMpi();
    ~IOMpi();
    int Cprintf();

private:
    const int BUFFER_SIZE = 10240;

    int Get_io_rank();
    MPI_Comm io_comm;
    int io_rank, my_rank, comm_size;
    int cprintf_tag;
    char *io_buff;
};


#endif //MPI3SNP_IOMpi_H
