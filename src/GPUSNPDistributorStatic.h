/*
 * GPUSNPDistributorStatic.h
 *
 *  Created on: 15/10/2014
 *      Author: jorge
 */

#ifndef GPUSNPDISTRIBUTORSTATIC_H_
#define GPUSNPDISTRIBUTORSTATIC_H_

#include "GPUSNPDistributor.h"

class GPUSNPDistributorStatic : public GPUSNPDistributor {
public:
        GPUSNPDistributorStatic(Options* options);
        ~GPUSNPDistributorStatic();

        inline void loadSNPSet(){
        	GPUSNPDistributor::loadSNPSet();
        	_initializeArrays();
        }

    	// Returns the number of pairs that will be computed (all the 3-way combinations)
    	// Also return the ids of the pairs
    	uint32_t getPairsSNPs(uint2* ids, uint64_t& totalAnal, uint16_t gpuId);

private:
        void _initializeArrays();

        uint32_t *_iterDoubleSnp1Arr;
        uint32_t *_iterDoubleSnp2Arr;
        bool *_morePairsArr;
};

#endif // GPUSNPDISTRIBUTORSTATIC_H_

