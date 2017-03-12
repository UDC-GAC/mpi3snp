//
// Created by christian on 10/03/17.
//

#include "IOMpi.h"

IOMpi::IOMpi() {
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    MPI_Comm_size(io_comm, &comm_size);
    MPI_Comm_rank(io_comm, &my_rank);
    io_rank = MPI_PROC_NULL;
    cprintf_tag = 0;
    pthread_mutex_init(&cprintf_mutex, NULL);
    io_buff = NULL;
}

IOMpi::~IOMpi() {
    pthread_mutex_destroy(&cprintf_mutex);
    if (io_buff != NULL){
        delete io_buff;
    }
}

int IOMpi::Get_io_rank() {
    // Determine IO process
    int mpi_io, flag, io_rank;
    MPI_Attr_get(MPI_COMM_WORLD, MPI_IO, &mpi_io, &flag);
    if (!flag) {
        // Attribute not cached
        io_rank = 0;
    } else if (mpi_io == MPI_PROC_NULL) {
        // No process can carry IO
        io_rank = MPI_PROC_NULL;
    } else if (mpi_io == MPI_ANY_SOURCE) {
        io_rank = 0;
    } else {
        // Multiple IO processes
        MPI_Allreduce(&mpi_io, &io_rank, 1, MPI_INT, MPI_MIN, io_comm);
    }
    return io_rank;
}

// Collective printf
int IOMpi::Cprintf(char *format, ...) {
    va_list args;
    int tag, i;

    if (io_rank == MPI_PROC_NULL){
        io_rank = Get_io_rank();
    }

    if (io_buff == NULL){
        io_buff = new char[BUFFER_SIZE];
    }

    pthread_mutex_lock(&cprintf_mutex);
    tag = cprintf_tag++;
    pthread_mutex_unlock(&cprintf_mutex);

    if (my_rank == io_rank){
        for (i=0; i<my_rank; i++){
            if (MPI_Recv(io_buff, BUFFER_SIZE, MPI_CHAR, i, tag, io_comm, NULL) != MPI_SUCCESS){
                return -1;
            }
            printf("Process %d > %s", i, io_buff);
            fflush(stdout);
        }
        va_start(args, format);
        vsprintf(io_buff, format, args);
        va_end(args);
        printf("Process %d > %s", my_rank, io_buff);
        for (i=my_rank+1; i<comm_size; i++){
            if (MPI_Recv(io_buff, BUFFER_SIZE, MPI_CHAR, i, tag, io_comm, NULL) != MPI_SUCCESS){
                return -1;
            }
            printf("Process %d > %s", i, io_buff);
            fflush(stdout);
        }
    } else {
        va_start(args, format);
        vsprintf(io_buff, format, args);
        va_end(args);
        if (MPI_Send(io_buff, strlen(io_buff) + 1, MPI_CHAR, io_rank, tag, io_comm) != MPI_SUCCESS){
            return -1;
        }
    }
    /* TODO:
     - return number of characters printed
    */
    return 0;
}