/*
 * SNPDistributor.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#ifndef SNPDISTRIBUTOR_H_
#define SNPDISTRIBUTOR_H_

#include "SNP.h"
#include "Options.h"
#include "MyFile.h"
#include "LineReader.h"
#include "MutualInfo.h"
#include "BoolVector.h"

class SNPDistributor {
public:
    SNPDistributor(Options *options);

    virtual ~SNPDistributor();

    // Loads all the SNPs in the distributor
    void loadSNPSet();

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
        return _numSnp;
    }

    inline uint16_t getNumCases() {
        return bv.falseCount();
    }

    inline uint16_t getNumCtrls() {
        return bv.trueCount();
    }

    inline vector<SNP *> getSnpSet() {
        return _snpSet;
    }

private:
    // Load the information of which individuals are control and case from the tfam file
    void _loadIndsClass();

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
    vector<SNP *> _snpSet;
    // File handler
    MyFilePt _fpTfam;
    MyFilePt _fpTped;
    MyFilePt _fpOut;
    LineReader *_lineReader;
    bool _withLock;

    // Variables to indicate if there are more analyses to perform
    bool _moreDouble;

    // Variables shared among all threads
    pthread_mutex_t _mutex;

    // Array to keep which values are cases
    // False cases, true controls
    BoolVector bv;

    // To store the dimensions of the problem
    uint32_t _numSnp;

    // Iterators for the SNPs
    uint32_t _iterDoubleSnp1;
    uint32_t _iterDoubleSnp2;
};

#endif /* SNPDISTRIBUTOR_H_ */
