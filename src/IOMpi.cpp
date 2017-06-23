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
        delete[] io_buff;
    }
}

int IOMpi::Get_io_rank() {
    // Determine IO process
    int mpi_io, flag, io_rank;

    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_IO, &mpi_io, &flag);
    if (!flag) {
        // Attribute not cached
        io_rank = DEFAULT_IO_PROC;
    } else if (mpi_io == MPI_PROC_NULL || mpi_io == MPI_ANY_SOURCE || comm_size < mpi_io) {
        io_rank = DEFAULT_IO_PROC;
    } else {
        // Multiple IO processes
        int err = MPI_Allreduce(&mpi_io, &io_rank, 1, MPI_INT, MPI_MIN, io_comm);
        if (err != MPI_SUCCESS){
            int error_class, len;
            char error_string[1024];
            MPI_Error_class(err, &error_class);
            MPI_Error_string(error_class, error_string, &len);
            std::cerr << error_string << std::endl;
            io_rank = DEFAULT_IO_PROC;
        }
    }
    return io_rank;
}

// Collective printf
int IOMpi::Cprintf(const char *format, ...) {
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

int IOMpi::Mprintf(const char *format, ...) {
    if (my_rank != 0){
        return 0;
    }

    if (io_buff == NULL){
        io_buff = new char[BUFFER_SIZE];
    }

    va_list args;
    va_start(args, format);
    vsprintf(io_buff, format, args);
    va_end(args);
    printf(io_buff);
    fflush(stdout);

    return 0;
}