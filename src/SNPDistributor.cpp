/*
 * SNPDistributor.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#include "SNPDistributor.h"

SNPDistributor::SNPDistributor(Options *options) {
    _options = options;
    _snpSet.reserve(DEFAULT_NUM_SNPS);
    bv.reserve(DEFAULT_NUM_INDS);

    _withLock = (_options->getNumThreads()) > 1;
    pthread_mutex_init(&_mutex, NULL);

    if ((_fpOut = myfopen(options->getOutFileName().c_str(), "wb")) == NULL) {
        Utils::exit("Out file: file %s could not be opened\n", options->getOutFileName().c_str());
    }

    if ((_fpTfam = myfopen(options->getTFAMFileName().c_str(), "r")) == NULL) {
        Utils::exit("TFAM file: file %s could not be opened\n", options->getTFAMFileName().c_str());
    }

    if ((_fpTped = myfopen(options->getTPEDFileName().c_str(), "r")) == NULL) {
        Utils::exit("TPED file: file %s could not be opened\n", options->getTPEDFileName().c_str());
    }

    _lineReader = new LineReader();

    _moreDouble = true;
}

SNPDistributor::~SNPDistributor() {
    if (dist != NULL){
        delete[] dist;
    }
    myfclose(_fpTfam);
    myfclose(_fpTped);
    myfclose(_fpOut);
}

void SNPDistributor::_loadIndsClass() {
    // Load the information from a TFAM file about the cases and controls
    int retValue;

    while ((retValue = _lineReader->readTFAMLine(_fpTfam, _snpSet.size())) >= 0) {
        bv.push_back(retValue);
    }

    //printf("Loaded information of %ld individuals (%ld/%ld cases/controls)\n", numInds, getNumCases(), getNumCtrls());
#ifdef DEBUG
    for(int i=0; i<numInds; i++){
        if(_indsClass[i]){
            Utils::log("control ");
        }
        else{
            Utils::log("case ");
        }
    }
    Utils::log("\n");
#endif
}

void SNPDistributor::loadSNPSet() {
    // Load the information of which individuals are cases and controls
    _loadIndsClass();

    SNP *readSNP = new(SNP);

    while (_lineReader->readTPEDLine(_fpTped, readSNP, _snpSet.size(), bv.size(), &bv[0])) {
        _snpSet.push_back(readSNP);
        readSNP = new(SNP);
    }
    _moreDouble = !_snpSet.empty();
    delete readSNP;

    dist_size = _snpSet.size() / _options->getNumProcesses() +
            (_options->getProcessId() < (_snpSet.size() % _options->getNumProcesses()));
    dist = new uint32_t[dist_size];
    dist[0] = _options->getProcessId();
    dist[1] = 2 * _options->getNumProcesses() - _options->getProcessId() - 1;
    for (int i = 2; i < dist_size; i++) {
        dist[i] = dist[i - 2] + 2 * _options->getNumProcesses();
    }

    dist_it = 0;
    _iterDoubleSnp1 = dist[dist_it++];
    _iterDoubleSnp2 = _iterDoubleSnp1 + 1;

#ifdef DEBUG
    int j;
    uint16_t casesAa, ctrlsAa;
    for(int i=0; i<_numSnp; i++){
        Utils::log("SNP %d:\n   name %s\n   cases 0", i, (char *)_snpSet[i]->_name);
        for(j=0; j<(_numCases/32+((_numCases%32)>0)); j++){
            Utils::log("%u ", _snpSet[i]->_case0Values[j]);
        }

        Utils::log("\n   cases 2 ");
        for(j=0; j<(_numCases/32+((_numCases%32)>0)); j++){
            Utils::log("%u ", _snpSet[i]->_case2Values[j]);
        }
        Utils::log("\n   controls 0 ");
        for(j=0; j<(_numCtrls/32+((_numCtrls%32)>0)); j++){
            Utils::log("%u ", _snpSet[i]->_ctrl0Values[j]);
        }

        Utils::log("\n   controls 2 ");
        for(j=0; j<(_numCtrls/32+((_numCtrls%32)>0)); j++){
            Utils::log("%u ", _snpSet[i]->_ctrl2Values[j]);
        }
    }
#endif
}

uint32_t SNPDistributor::_getPairsSNPsNoLock(uint32_t *ids) {
    uint16_t iter_block = 0;

    if (_moreDouble) {
        while (dist_it < dist_size) {
            while (_iterDoubleSnp2 < _snpSet.size() && iter_block < NUM_PAIRS_BLOCK) {
                ids[2 * iter_block] = _iterDoubleSnp1;
                ids[2 * iter_block + 1] = _iterDoubleSnp2++;
                iter_block++;
            }

            if (iter_block == NUM_PAIRS_BLOCK) {
                return iter_block;
            } else {
                _iterDoubleSnp1 = dist[dist_it++];
                _iterDoubleSnp2 = _iterDoubleSnp1 + 1;
            }
        }
        _moreDouble = false;
    }

    return iter_block;
}

uint32_t SNPDistributor::_getPairsSNPsLock(uint32_t *ids) {
    uint16_t iter_block = 0;

    _lock();
    iter_block = _getPairsSNPsNoLock(ids);
    _unlock();

    return iter_block;
}
