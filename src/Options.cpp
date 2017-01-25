/*
 * Options.cpp
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#include "Options.h"

Options::Options(int *argc, char ***argv) {
    _setDefaults();

    _c = argc;
    _v = argv;
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
            _numCPUs);

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
    _numCPUs = 1;
    _numOutputs = 10;
}

bool Options::parse() {
    int argc = *_c;
    char **argv = *_v;
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
                Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        } else if (!strcmp(argv[argind], "-rf")) {
            argind++;
            if (argind < argc) {
                _tfamFileName = argv[argind];
                argind++;
            } else {
                Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
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
                Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
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
                _numCPUs = intVal;
            } else {
                Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
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
                Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
                return false;
            }
        }
    }

    Utils::log("Number of CPU threads: %d\n", _numCPUs);
    Utils::log("\nNumber of outputs: %hu\n", _numOutputs);
    Utils::log("Number of pairs by block: %hu\n", NUM_PAIRS_BLOCK);

    if (!strcmp(_tfamFileName.c_str(), "") || !strcmp(_tpedFileName.c_str(), "")) {
        Utils::log("Input files not specified!!!\n");
        return false;
    }

    if (!strcmp(_outFileName.c_str(), "")) {
        Utils::log("Output file not specified!!!\n");
        return false;
    }

    return true;
}
