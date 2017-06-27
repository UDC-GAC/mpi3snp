/*
 * SNPDistributor.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#ifndef SNPDISTRIBUTOR_H_
#define SNPDISTRIBUTOR_H_

#include <mpi.h>
#include "SNP2.h"
#include "MyFile.h"
#include "LineReader.h"
#include "BoolVector.h"

class SNPDistributor {
public:
    SNPDistributor(int num_proc, int proc_id, int num_thread, std::string tped, std::string tfam);

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

    inline uint32_t getNumSnp() {
        return _snpSet.size();
    }

    inline uint16_t getNumCases() {
        return bv.falseCount();
    }

    inline uint16_t getNumCtrls() {
        return bv.trueCount();
    }

    inline vector<SNP2 *> &getSnpSet() {
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

    vector<SNP2 *> _snpSet;
    // File handler
    MyFilePt _fpTfam;
    MyFilePt _fpTped;
    LineReader *_lineReader;
    bool _withLock;

    // Variables to indicate if there are more analyses to perform
    bool _moreDouble;

    // Variables shared among all threads
    pthread_mutex_t _mutex;

    // Array to keep which values are cases
    // False cases, true controls
    BoolVector bv;

    // Iterators for the SNPs
    int dist_size, dist_it;
    uint32_t *dist = NULL;
    uint32_t _iterDoubleSnp1;
    uint32_t _iterDoubleSnp2;

    int num_proc, proc_id, num_threads;
};

#endif /* SNPDISTRIBUTOR_H_ */
