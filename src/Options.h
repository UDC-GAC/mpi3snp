/*
 * Options.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <mpi.h>
#include "Utils.h"
#include "IOMpi.h"

class Options {
public:
    Options();

    virtual ~Options();

    inline int getNumThreads() {
        return _numThreads;
    }

    inline int getNumProcesses() {
        return _numProcesses;
    }

    inline int getProcessId() {
        return _processId;
    }

    inline string &getTPEDFileName() {
        return _tpedFileName;
    }

    inline string &getTFAMFileName() {
        return _tfamFileName;
    }

    inline string &getOutFileName() {
        return _outFileName;
    }

    inline uint16_t getNumOutputs() {
        return _numOutputs;
    }

    // To parse the command line
    bool parse(int argc, char **argv);

    // To print the help for the users
    void printUsage();

private:
    string _tpedFileName;
    string _tfamFileName;
    string _outFileName;
    int _numThreads;
    int _numProcesses;
    int _processId;
    uint16_t _numOutputs;

    /*private member functions*/
    void _setDefaults();
};

#endif /* OPTIONS_H_ */
