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
#include "Block.h"
#include "IOMpi.h"

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
        return _snpSet.size();
    }

    inline uint16_t getNumCases() {
        return _indsClass.falseCount();
    }

    inline uint16_t getNumCtrls() {
        return _indsClass.trueCount();
    }

    inline vector<SNP *> getSnpSet() {
        return _snpSet;
    }

    inline BoolVector getClassSet() {
        return _indsClass;
    }

    inline void setSNPBlocks(std::vector<std::vector<Block *>> blockList) {
        this->blockList = blockList;
        blockIt = this->blockList.begin();
        index1 = (*blockIt)[0]->x;
        index1Lim = index1 + (*blockIt)[0]->xlen;
        index2 = (*blockIt)[1]->x;
        index2Lim = index2 + (*blockIt)[1]->xlen;

        isDiagonal = index1 == index2;
        if (isDiagonal) {
            index2 = index1 + 1;
        }
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
    std::vector<SNP *> _snpSet;
    std::vector<std::vector<Block *>> blockList;
    // Iterators for the SNPs
    size_t index1, index1Lim, index2, index2Lim;
    std::vector<std::vector<Block *>>::iterator blockIt;
    bool isDiagonal;

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
    BoolVector _indsClass;
};

#endif /* SNPDISTRIBUTOR_H_ */
