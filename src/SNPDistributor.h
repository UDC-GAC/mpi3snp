/*
 * SNPDistributor.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 *
 *  Modified on: December 29, 2016
 *           by: Christian Ponte Fernández
 */

#ifndef SNPDISTRIBUTOR_H_
#define SNPDISTRIBUTOR_H_

#include "SNP.h"
#include "Options.h"
#include "MyFile.h"
#include "LineReader.h"
#include "MutualInfo.h"

class SNPDistributor {
public:
    // TODO: using struct with vector and counters instead of creating a specific class
    static typedef struct {
        vector<char> classVector;
        uint16_t numCases;
        uint16_t numCtrls;
    } ClassSet_t;

    SNPDistributor(Options *options);

    SNPDistributor(Options *options, vector<SNP *> snpSet, ClassSet_t classSet, uint32_t lowerLim, uint32_t upperLim);

    virtual ~SNPDistributor();

    // Returns the number of pairs that will be computed (all the 3-way combinations)
    // Also return the ids of the pairs
    inline uint32_t getPairsSNPs(uint32_t *ids) {
        if (_withLock) {
            return _getPairsSNPsLock(ids);
        } else {
            return _getPairsSNPsNoLock(ids);
        }
    }

    inline void printMI(MutualInfo *info, uint16_t numOutputs) {
        for (int i = numOutputs - 1; i >= 0; i--) {
            MutualInfo auxInfo = info[i];
            _lineReader->printMI(_fpOut, auxInfo._id1, auxInfo._id2, auxInfo._id3, auxInfo._mutualInfoValue);
        }
    }

    inline uint32_t getNumSnp() {
        return _snpSet.size();
    }

    inline uint16_t getNumCases() {
        return _classSet.numCases;
    }

    inline uint16_t getNumCtrls() {
        return _classSet.numCtrls;
    }

    inline vector<SNP *> getSnpSet() {
        return _snpSet;
    }

    inline ClassSet_t getClassSet(){
        return _classSet;
    }

private:
    // Loads all the SNPs in the distributor
    void _loadSNPSet(MyFilePt fpTped);

    // Load the information of which individuals are control and case from the tfam file
    void _loadIndsClass(MyFilePt fpTfam);

    // Functions to work with mutex
    inline void _lock() {
        if (_withLock) {
            pthread_mutex_lock(&_mutex);
        }
    }

    inline void _unlock() {
        if (_withLock) {
            pthread_mutex_unlock(&_mutex);
        }
    }

    uint32_t _getPairsSNPsLock(uint32_t *ids);

    uint32_t _getPairsSNPsNoLock(uint32_t *ids);

    Options *_options;
    // File handler
    LineReader *_lineReader;
    MyFilePt _fpOut;

    // Variables shared among all threads
    pthread_mutex_t _mutex;
    bool _withLock;
    vector<SNP *> _snpSet;
    uint32_t _upperLim;
    // Variables to indicate if there are more analyses to perform
    bool _moreDouble;

    // Array to keep which values are cases
    // False cases, true controls
    ClassSet_t _classSet;

    // Iterators for the SNPs
    uint32_t _iterDoubleSnp1;
    uint32_t _iterDoubleSnp2;
};

#endif /* SNPDISTRIBUTOR_H_ */
