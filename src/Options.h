/*
 * Options.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "Utils.h"

class Options {
public:
    Options(int *argc, char ***argv);

    virtual ~Options();

    inline int getNumCPUs() {
        return _numCPUs;
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

    inline int *get_argc(){
        return _c;
    }

    inline char ***get_argv(){
        return _v;
    }

    // To parse the command line
    bool parse();

    // To print the help for the users
    void printUsage();

private:
    int *_c;
    char ***_v;

    string _tpedFileName;
    string _tfamFileName;
    string _outFileName;
    int _numCPUs;
    uint16_t _numOutputs;

    /*private member functions*/
    void _setDefaults();
};

#endif /* OPTIONS_H_ */
