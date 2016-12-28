/*
 * SNPDistributor.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 *
 *  Modified on: December 28, 2016
 *           by: Christian Ponte Fernández
 */

#include "SNPDistributor.h"

SNPDistributor::SNPDistributor(Options *options) {
    _options = options;

    _withLock = (_options->getNumCPUs()) > 1;
    pthread_mutex_init(&_mutex, NULL);

    _lineReader = new LineReader();
    MyFilePt fpTfam, fpTped;

    if ((_fpOut = myfopen(options->getOutFileName().c_str(), "wb")) == NULL) {
        Utils::exit("Out file: file %s could not be opened\n", options->getOutFileName().c_str());
    }

    if ((fpTfam = myfopen(options->getTFAMFileName().c_str(), "r")) == NULL) {
        Utils::exit("TFAM file: file %s could not be opened\n", options->getTFAMFileName().c_str());
    }

    if ((fpTped = myfopen(options->getTPEDFileName().c_str(), "r")) == NULL) {
        Utils::exit("TPED file: file %s could not be opened\n", options->getTPEDFileName().c_str());
    }

    _classSet.classVector.resize(DEFAULT_NUM_INDS);
    _classSet.numCases = 0;
    _classSet.numCtrls = 0;
    _loadIndsClass(fpTfam);
    myfclose(fpTfam);

    _snpSet.resize(DEFAULT_NUM_SNPS);
    _loadSNPSet(fpTped);
    myfclose(fpTped);

    _upperLim = _snpSet.size();
    _iterDoubleSnp1 = 0;
    _iterDoubleSnp2 = _iterDoubleSnp1 + 1;
}

SNPDistributor::SNPDistributor(Options *options, vector<SNP *> snpSet, ClassSet_t classSet, uint32_t lowerLim,
                               uint32_t upperLim) {
    _options = options;

    _withLock = (_options->getNumCPUs()) > 1;
    pthread_mutex_init(&_mutex, NULL);

    _lineReader = new LineReader();

    if ((_fpOut = myfopen(options->getOutFileName().c_str(), "wb")) == NULL) {
        Utils::exit("Out file: file %s could not be opened\n", options->getOutFileName().c_str());
    }

    _classSet = classSet;

    _snpSet = snpSet;

    _upperLim = upperLim;
    _iterDoubleSnp1 = lowerLim;
    _iterDoubleSnp2 = _iterDoubleSnp1 + 1;
}

SNPDistributor::~SNPDistributor() {
    myfclose(_fpOut);
}

void SNPDistributor::_loadSNPSet(MyFilePt fpTped) {
    SNP *readSNP = new(SNP);

    // TODO: LineReader requires bool * instead of vector<uint8_t> | vector<char>
    while (_lineReader->readTPEDLine(fpTped, readSNP, _snpSet.size(), _classSet.classVector.size(),
                                     (bool *) &_classSet.classVector[0])) {
        _snpSet.push_back(readSNP);

        readSNP = new(SNP);
    }

    delete readSNP;

    _moreDouble = _snpSet.size();
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

void SNPDistributor::_loadIndsClass(MyFilePt fpTfam) {
    // Load the information from a TFAM file about the cases and controls
    int retValue;
    _classSet.numCases = 0;
    _classSet.numCtrls = 0;

    while ((retValue = _lineReader->readTFAMLine(fpTfam, _classSet.classVector.size())) >= 0) {

        if (retValue) {
            _classSet.numCtrls++;
            _classSet.classVector.push_back(true);
        } else {
            _classSet.numCases++;
            _classSet.classVector.push_back(false);
        }
    }

    Utils::log("Loaded information of %ld individuals (%ld/%ld cases/controls)\n", _classSet.classVector.size(),
               _classSet.numCases, _classSet.numCtrls);
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


uint32_t SNPDistributor::_getPairsSNPsNoLock(uint32_t *ids) {
    uint32_t id1;
    uint32_t id2;

    if (_moreDouble) {
        uint16_t iter_block = 0;
        while (iter_block < NUM_PAIRS_BLOCK) {

            id1 = _iterDoubleSnp1;
            id2 = _iterDoubleSnp2;

            ids[2 * iter_block] = id1;
            ids[2 * iter_block + 1] = id2;

            iter_block++;

            // Look for the next pair
            if (_iterDoubleSnp2 == _snpSet.size() - 1) {
                _iterDoubleSnp1++;
                _iterDoubleSnp2 = _iterDoubleSnp1 + 1;
            } else {
                _iterDoubleSnp2++;
            }

            if (_iterDoubleSnp1 == _upperLim - 1) { // We have finished to compute the block
                _moreDouble = false;
                break;
            }
        }

        return iter_block;
    }
    return 0;
}

uint32_t SNPDistributor::_getPairsSNPsLock(uint32_t *ids) {
    uint16_t iter_block;

    _lock();
    iter_block = _getPairsSNPsNoLock(ids);
    _unlock();

    return iter_block;
}
