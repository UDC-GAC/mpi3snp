//
// Created by christian on 10/03/17.
//

#include "IOMpi.h"

IOMpi::IOMpi() {
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    MPI_Comm_size(io_comm, &comm_size);
    MPI_Comm_rank(io_comm, &my_rank);
    io_rank = Get_io_rank();
    cprintf_tag = 0;
    pthread_mutex_init(&cprintf_mutex, NULL);
}

IOMpi::~IOMpi() {
    pthread_mutex_destroy(&cprintf_mutex);
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
        if (err != MPI_SUCCESS) {
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


int IOMpi::Cfprintf_list(std::ostream &ostream, const char *format, va_list &list){
    int tag, charcount = 0;

    pthread_mutex_lock(&cprintf_mutex);
    tag = cprintf_tag++;
    pthread_mutex_unlock(&cprintf_mutex);

    if (my_rank == io_rank) {
        MPI_Message message;
        MPI_Status status;
        int count;
        for (int i = 0; i < comm_size; i++) {
            if (i == my_rank) {
                const int str_size = snprintf(nullptr, 0, format, list);
                char str[str_size + 1];
                vsprintf(str, format, list);
                ostream << "Process " + std::to_string(my_rank) + " > " << str;
                charcount += str_size;
            } else {
                if (MPI_Mprobe(i, tag, io_comm, &message, &status) != MPI_SUCCESS) {
                    return -1;
                }
                MPI_Get_count(&status, MPI_CHAR, &count);
                char str[count];
                if (MPI_Mrecv(str, count, MPI_CHAR, &message, NULL) != MPI_SUCCESS) {
                    return -1;
                }
                ostream << "Process " + std::to_string(i) + " > " << str;
                charcount += count;
            }
        }
        std::flush(ostream);
    } else {
        const int str_size = snprintf(nullptr, 0, format, list);
        char str[str_size + 1];
        vsprintf(str, format, list);
        if (MPI_Send(str, str_size + 1, MPI_CHAR, io_rank, tag, io_comm) != MPI_SUCCESS) {
            return -1;
        }
        charcount = str_size;
    }
    return charcount;
}

int IOMpi::Mfprintf_list(std::ostream &ostream, const char *format, va_list &list) {
    if (my_rank != io_rank) {
        return 0;
    }

    const int str_size = snprintf(nullptr, 0, format, list);
    char str[str_size + 1];
    vsprintf(str, format, list);
    ostream << str;
    std::flush(ostream);
    return str_size;
}