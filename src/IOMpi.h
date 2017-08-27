//
// Created by christian on 10/03/17.
//

#ifndef MPI3SNP_IOMpi_H
#define MPI3SNP_IOMpi_H

#include <cstdarg>
#include <iostream>
#include <pthread.h>

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
    enum Level {
        E = 0, // Error
        N = 1, // None
        B = 2, // Benchmarking
        D = 3  // Debug
    };

    static IOMpi &Instance() {
        static IOMpi inst;
        return inst;
    }

    static void Set_print_level(Level l) {
        Instance().level = l;
    }

    // delete copy and move constructors and assign operators
    IOMpi(IOMpi const &) = delete;             // Copy construct
    IOMpi(IOMpi &&) = delete;                  // Move construct
    IOMpi &operator=(IOMpi const &) = delete;  // Copy assign
    IOMpi &operator=(IOMpi &&) = delete;      // Move assign

    template<Level l>
    inline int Cprintf(const char *format, ...) {
        int c = 0;
        if (l <= level) {
            va_list list;
            va_start(list, format);
            c = Cfprintf_list(std::cout, format, list);
            va_end(list);
        }
        return c;
    }

    template<Level l>
    inline int Cfprintf(std::ostream &ostream, const char *format, ...) {
        int c = 0;
        if (l <= level) {
            va_list list;
            va_start(list, format);
            c = Cfprintf_list(ostream, format, list);
            va_end(list);
        }
        return c;
    }

    template<Level l>
    inline int Mprintf(const char *format, ...) {
        int c = 0;
        if (l <= level) {
            va_list list;
            va_start(list, format);
            c = Mfprintf_list(std::cout, format, list);
            va_end(list);
        }
        return c;
    }

    template<Level l>
    inline int Mfprintf(std::ostream &ostream, const char *format, ...) {
        int c = 0;
        if (l <= level) {
            va_list list;
            va_start(list, format);
            c = Mfprintf_list(ostream, format, list);
            va_end(list);
        }
        return c;
    }

protected:
    /* Methods */
    IOMpi();

    ~IOMpi();

    int Get_io_rank();

    int Cfprintf_list(std::ostream &ostream, const char *format, va_list &list);

    int Mfprintf_list(std::ostream &ostream, const char *format, va_list &list);

    /* Attributes */
    static const int DEFAULT_IO_PROC = 0;

    void * void_io_comm;
    int io_rank, my_rank, comm_size;
    int cprintf_tag;
    pthread_mutex_t cprintf_mutex;
    Level level;
};

#endif //MPI3SNP_IOMpi_H