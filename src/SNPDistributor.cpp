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

    _withLock = (_options->getNumCPUs()) > 1;
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

    _indsClass.reserve(DEFAULT_NUM_INDS);
    _lineReader = new LineReader();

    _moreDouble = true;

    index1 = 0;
    index2 = 1;
}

SNPDistributor::~SNPDistributor() {

    myfclose(_fpTfam);
    myfclose(_fpTped);
    myfclose(_fpOut);
}

void SNPDistributor::_loadIndsClass() {
    // Load the information from a TFAM file about the cases and controls
    int retValue;

    while ((retValue = _lineReader->readTFAMLine(_fpTfam, _indsClass.size())) >= 0) {
        _indsClass.push_back(retValue);
    }

    IOMpi::Instance().Cprintf("Loaded information of %ld individuals (%ld/%ld cases/controls)\n", _indsClass.size(),
               _indsClass.falseCount(), _indsClass.trueCount());
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

    SNP *readSNP;

    while (1) {
        readSNP = new(SNP);

        if (!_lineReader->readTPEDLine(_fpTped, readSNP, _snpSet.size(), _indsClass.size(), &_indsClass[0])) {
            delete readSNP;
            break;
        }

        _snpSet.push_back(readSNP);
    }

    if (!_snpSet.size()) {
        _moreDouble = false;
    }

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
    if (_moreDouble) {
        uint16_t iter_block = 0;
        while (iter_block < NUM_PAIRS_BLOCK) {
            ids[2 * iter_block] = index1;
            ids[2 * iter_block + 1] = index2;

            iter_block++;
            // Look for the next pair
            if (++index2 == index2Lim) {
                // Update index1 and restart index2
                index1++;
                if (isDiagonal) {
                    if (index1 == index1Lim - 1){
                        index1++;
                    } else {
                        index2 = index1 + 1;
                    }
                } else {
                    index2 = (*blockIt)[1]->x;
                }

                // Update blockIt based on index1
                if (index1 == index1Lim){
                    if (++blockIt != blockList.end()){
                        index1 = (*blockIt)[0]->x;
                        index1Lim = index1 + (*blockIt)[0]->xlen;
                        index2 = (*blockIt)[1]->x;
                        index2Lim = index2 + (*blockIt)[1]->xlen;
                    } else {
                        _moreDouble = false;
                        break;
                    }
                }
            }
        }
        return iter_block;
    }
    return 0;
}

uint32_t SNPDistributor::_getPairsSNPsLock(uint32_t *ids) {
    uint16_t iter_block = 0;

    _lock();
    iter_block = _getPairsSNPsNoLock(ids);
    _unlock();

    return iter_block;
}
