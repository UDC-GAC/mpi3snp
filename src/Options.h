/*
 * Options.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <cstdint>
#include <vector>

class Options {
public:
    Options();

    Options(int argc, char *argv[]);

    virtual ~Options();

    inline std::vector<uint16_t> Get_GPU_Ids() {
        return _GPUIds;
    }

    inline int getNumProcesses() {
        return _numProcesses;
    }

    inline int getProcessId() {
        return _processId;
    }

    inline char *&getTPEDFileName() {
        return _tpedFileName;
    }

    inline char *&getTFAMFileName() {
        return _tfamFileName;
    }

    inline char *&getOutFileName() {
        return _outFileName;
    }

    inline uint16_t getNumOutputs() {
        return _numOutputs;
    }

    inline bool isHeteroGPUs() {
        return _heteroGPUs;
    }

    inline bool isMI() {
        return _filterMI;
    }

    // To parse the command line
    bool parse(int argc, char *argv[]);

    // To print the help for the users
    void printUsage();

private:
    char *_tpedFileName;
    char *_tfamFileName;
    char *_outFileName;
    int _numProcesses;
    int _processId;
    std::vector<uint16_t> _GPUIds;
    uint16_t _numOutputs;
    bool _heteroGPUs;
    bool _filterMI;

    /*private member functions*/
    void _setDefaults();
};

#endif /* OPTIONS_H_ */
