/*
 * Options.cpp
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#include "Options.h"

Options::Options() {
    _setDefaults();

    MPI_Comm_size(MPI_COMM_WORLD, &_numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &_processId);
}

Options::~Options() {

}

void Options::printUsage() {
    fprintf(stderr,
            "--------------------------------------------------------------------------------\n\n");
    fprintf(stderr,
            "3SNPInt (v%s) is a tool to detect 3-way interaction in Genome Wide Association datasets using"
                    "multicore CPUs\n", VERSION);

    fprintf(stderr,
            "--------------------------------------------------------------------------------\n\n");

    fprintf(stderr, "Usage: 3SNPInt-Exh -rp srcTPEDFile -rf srcTFAMFile -o outFile [options]\n");

    // Input
    fprintf(stderr, "\tsrcTPEDFile: the file name for the SNPs inputs with -tped format\n");
    fprintf(stderr, "\tsrcTFAMFile: the file name for the individuals information with -tfam format\n");

    // Output
    fprintf(stderr, "Output:\n");
    fprintf(stderr, "\toutFile: output file name\n");

    // Compute
    fprintf(stderr, "Computational options:\n");
    fprintf(stderr, "\t-t <int> (CPU threads for computation, default = %d)\n",
            _numThreads);

    fprintf(stderr, "\t-no <int> (number of outputs, default = %hu)\n",
            _numOutputs);

    fprintf(stderr, "Others:\n");
    fprintf(stderr, "\t-h <print out the usage of the program)\n");
}

void Options::_setDefaults() {
    // Empty string means outputing to STDOUT
    _tpedFileName = "";
    _tfamFileName = "";
    _outFileName = "";
    _numThreads = 1;
    _numOutputs = 10;
}

bool Options::parse(int argc, char **argv) {
    int intVal, i, j;
    uint16_t uintVal;
    int argind = 1;

    if (argc < 2) {
        return false;
    }

    // Check the help
    if (!strcmp(argv[argind], "-h") || !strcmp(argv[argind], "-?")) {
        return false;
    }

    // For the other options
    while (argind < argc) {
        // Input file
        if (!strcmp(argv[argind], "-rp")) {
            argind++;
            if (argind < argc) {
                _tpedFileName = argv[argind];
                argind++;
            } else {
                IOMpi::Instance().Mprintf("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        } else if (!strcmp(argv[argind], "-rf")) {
            argind++;
            if (argind < argc) {
                _tfamFileName = argv[argind];
                argind++;
            } else {
                IOMpi::Instance().Mprintf("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        } else if (!strcmp(argv[argind], "-o")) {
            argind++;
            if (argind < argc) {
                if (strlen(argv[argind]) > 0) {
                    _outFileName = argv[argind];
                    argind++;
                }
            } else {
                IOMpi::Instance().Mprintf("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        } else if (!strcmp(argv[argind], "-t")) {
            intVal = 1;
            argind++;
            if (argind < argc) {
                sscanf(argv[argind], "%d", &intVal);
                if (intVal < 1)
                    intVal = 1;

                argind++;
                _numThreads = intVal;
            } else {
                IOMpi::Instance().Mprintf("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        } else if (!strcmp(argv[argind], "-no")) {
            uintVal = 1;
            argind++;
            if (argind < argc) {
                sscanf(argv[argind], "%hu", &uintVal);
                argind++;
                _numOutputs = uintVal;
            } else {
                IOMpi::Instance().Mprintf("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        }
    }

    IOMpi::Instance().Mprintf("Number of MPI processes: %d\n", _numProcesses);
    IOMpi::Instance().Mprintf("Number of CPU threads: %d\n", _numThreads);
    IOMpi::Instance().Mprintf("Number of outputs: %hu\n", _numOutputs);
    IOMpi::Instance().Mprintf("Number of pairs by block: %hu\n", NUM_PAIRS_BLOCK);

    if (!strcmp(_tfamFileName.c_str(), "") || !strcmp(_tpedFileName.c_str(), "")) {
        IOMpi::Instance().Mprintf("Input files not specified!!!\n");
        return false;
    }

    if (!strcmp(_outFileName.c_str(), "")) {
        IOMpi::Instance().Mprintf("Output file not specified!!!\n");
        return false;
    }

    return true;
}
