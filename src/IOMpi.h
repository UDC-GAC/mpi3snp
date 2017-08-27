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
 *  singleton initialization (C++11 or superior). The object destruction is handled at the end of the program.
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
    inline int print(const std::string &s) {
        if (l <= level) {
            sprint_nolr(std::cout, s);
        }
    }

    template<Level l>
    inline int mprint(const std::string &s) {
        if (l <= level && my_rank == io_rank) {
            sprint_nolr(std::cout, s);
        }
    }

    template<Level l>
    inline int cprint(const std::string &s) {
        if (l <= level) {
            scprint_nol(std::cout, s);
        }
    }

    template<Level l>
    inline void sprint(std::ostream &ostream, const std::string &s) {
        if (l <= level) {
            sprint_nolr(ostream, s);
        }
    }

    template<Level l>
    inline int smprint(std::ostream &ostream, const std::string &s) {
        if (l <= level && my_rank == io_rank) {
            sprint_nolr(ostream, s);
        }
    }

    template<Level l>
    inline void scprint(std::ostream &ostream, const std::string &s) {
        if (l <= level) {
            scprint_nol(ostream, s);
        }
    }

protected:
    IOMpi();

    ~IOMpi();

    int Get_io_rank();

    void scprint_nol(std::ostream &ostream, const std::string &s);

    void sprint_nolr(std::ostream &ostream, const std::string &s);

    static const int DEFAULT_IO_PROC = 0;

    void * void_io_comm;
    int io_rank, my_rank, comm_size;
    int cprintf_tag;
    pthread_mutex_t cprintf_mutex;
    Level level;
};

#endif //MPI3SNP_IOMpi_H