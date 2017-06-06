/*
 * Options.cpp
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#include "Options.h"
#include "Macros.h"
#include "IOMpi.h"

Options::Options() {
    _setDefaults();

    MPI_Comm_size(MPI_COMM_WORLD, &_numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &_processId);
}

Options::Options(int argc, char *argv[]) {
    _setDefaults();

    MPI_Comm_size(MPI_COMM_WORLD, &_numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &_processId);

    // Parse the command line
    parse(argc, argv);
}

Options::~Options() {

    if (_GPUIds) {
        delete[] _GPUIds;
    }

}

void Options::printUsage() {
    fprintf(stderr,
            "--------------------------------------------------------------------------------\n\n");
    fprintf(stderr,
            "GPU3SNP (v%s) is a tool to detect 3-way interaction in Genome Wide Association datasets using"
                    "GPUs\n", VERSION);

    fprintf(stderr,
            "--------------------------------------------------------------------------------\n\n");

    fprintf(stderr, "Usage: GPU3SNP -rp srcTPEDFile -rf srcTFAMFile -o outFile [options]\n");

    // Input
    fprintf(stderr, "\tsrcTPEDFile: the file name for the SNPs inputs with -tped format\n");
    fprintf(stderr, "\tsrcTFAMFile: the file name for the individuals information with -tfam format\n");

    // Output
    fprintf(stderr, "Output:\n");
    fprintf(stderr, "\toutFile: output file name\n");

    // Compute
    fprintf(stderr, "Computational options:\n");
    fprintf(stderr, "\t-t <int> (number of GPUs for computation, default = %d)\n",
            _numGPUs);

    fprintf(stderr,
            "\t-g <int> [int] ... (index of the GPUs. There must be the same number of values as GPUs indicated with -t. Default = %d)\n",
            _GPUIds[0]);

    fprintf(stderr, "\t-no <int> (number of outputs, default = %hu)\n",
            _numOutputs);

    fprintf(stderr,
            "\t-het (to indicate that there are GPUs of different types and a dynamic distribution is applied)\n");

    fprintf(stderr, "\t-ig (apply Information Gain for filtering)\n");

    fprintf(stderr, "\t-mi (apply Mutual Information for filtering, by default)\n");

    fprintf(stderr, "Others:\n");
    fprintf(stderr, "\t-h <print out the usage of the program)\n");
}

void Options::_setDefaults() {
    // Empty string means outputing to STDOUT
    _tpedFileName = "";
    _tfamFileName = "";
    _outFileName = "";
    _numGPUs = 1;

    _GPUIds = new uint16_t[1];
    _GPUIds[0] = 0;

    _numOutputs = 10;
    _heteroGPUs = false;
    _filterMI = true;
}

bool Options::parse(int argc, char *argv[]) {
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

    // Print out the command line
    fprintf(stderr, "Command: ");
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, "%s ", argv[i]);
    }
    fputc('\n', stderr);

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
                _numGPUs = intVal;
                delete[] _GPUIds;
                _GPUIds = new uint16_t[_numGPUs];
            } else {
                IOMpi::Instance().Mprintf("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        } else if (!strcmp(argv[argind], "-g")) {
            argind++;
            for (i = 0; i < _numGPUs; i++) {
                if (argind < argc) {
                    sscanf(argv[argind], "%d", &intVal);
                    for (j = 0; j < i; j++) {
                        if (_GPUIds[j] == intVal) {
                            IOMpi::Instance().Mprintf("repeating GPU indexes not allowed\n");
                            return false;
                        }
                    }
                    _GPUIds[i] = intVal;
                } else {
                    IOMpi::Instance().Mprintf("not specify enough values for the parameter -g\n");
                    return false;
                }
                argind++;
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
        } else if (!strcmp(argv[argind], "-het")) {
            _heteroGPUs = true;
            argind++;
        } else if (!strcmp(argv[argind], "-ig")) {
            _filterMI = false;
            argind++;
        } else if (!strcmp(argv[argind], "-mi")) {
            _filterMI = true;
            argind++;
        }
    }

    if (_filterMI) {
        IOMpi::Instance().Mprintf("Applying Mutual Information\n");
    } else {
        IOMpi::Instance().Mprintf("Applying Information Gain\n");
    }
    IOMpi::Instance().Mprintf("Number of GPUs: %d\n", _numGPUs);
    IOMpi::Instance().Mprintf("GPU index:");
    for (i = 0; i < _numGPUs; i++) {
        IOMpi::Instance().Mprintf(" %d", _GPUIds[i]);
    }
    IOMpi::Instance().Mprintf("\nNumber of outputs: %hu\n", _numOutputs);
    IOMpi::Instance().Mprintf("Number of pairs by block: %lu\n", NUM_PAIRS_BLOCK);
    if (_numGPUs > 1) {
        if (_heteroGPUs) {
            IOMpi::Instance().Mprintf("Dynamic distribution among GPUs\n");
        } else {
            IOMpi::Instance().Mprintf("Static distribution among GPUs\n");
        }
    }

    if (!strcmp(_tfamFileName, "") || !strcmp(_tpedFileName, "")) {
        IOMpi::Instance().Mprintf("Input files not specified!!!\n");
        return false;
    }

    if (!strcmp(_outFileName, "")) {
        IOMpi::Instance().Mprintf("Output file not specified!!!\n");
        return false;
    }

    return true;
}
