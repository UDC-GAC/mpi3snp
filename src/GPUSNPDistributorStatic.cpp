/*
 * GPUSNPDistributorStatic.cpp
 *
 *  Created on: 15/10/2014
 *      Author: jorge
 */

#include "GPUSNPDistributorStatic.h"

GPUSNPDistributorStatic::GPUSNPDistributorStatic(Options* options):GPUSNPDistributor(options) {

	_iterDoubleSnp1Arr = new uint32_t[options->getNumGPUs()];
	_iterDoubleSnp2Arr = new uint32_t[options->getNumGPUs()];
	_morePairsArr = new bool[options->getNumGPUs()];

	for(int i=0; i<options->getNumGPUs(); i++){
		_morePairsArr[i] = true;
	}
}

GPUSNPDistributorStatic::~GPUSNPDistributorStatic(){

	if(_iterDoubleSnp1Arr){
		delete[] _iterDoubleSnp1Arr;
	}
	if(_iterDoubleSnp2Arr){
		delete[] _iterDoubleSnp2Arr;
	}
	if(_morePairsArr){
		delete[] _morePairsArr;
	}
}

void GPUSNPDistributorStatic::_initializeArrays(){

	uint32_t id1 = 0;
	uint32_t id2 = 1;

	for(int i=0; i<_options->getNumGPUs(); i++){

		_iterDoubleSnp1Arr[i] = id1;
		_iterDoubleSnp2Arr[i] = id2;

		// Look for the next pair
		if(id2 == _numSnp-2){
			id1++;
			id2 = id1+1;
		}
		else{
			id2++;
		}
	}
}

uint32_t GPUSNPDistributorStatic::getPairsSNPs(uint2* ids, uint64_t& totalAnal, uint16_t gpuId){

	uint32_t id1;
	uint32_t id2;

	// It applies a cyclic distribution among the pairs
	if(_morePairsArr[gpuId]){
		uint16_t iter_block = 0;
		while(iter_block < NUM_PAIRS_BLOCK){

			id1 = _iterDoubleSnp1Arr[gpuId];
			id2 = _iterDoubleSnp2Arr[gpuId];

			totalAnal += 1;//_numSnp-id2-1;

			ids[iter_block].x = id1;
			ids[iter_block].y = id2;

			iter_block++;

			// Look for the next pair
			if(id2 >= _numSnp-1-_options->getNumGPUs()){
				_iterDoubleSnp1Arr[gpuId]++;
				_iterDoubleSnp2Arr[gpuId] = _iterDoubleSnp1Arr[gpuId]+1+gpuId;
			}
			else{
				_iterDoubleSnp2Arr[gpuId] += _options->getNumGPUs();
			}

			if(_iterDoubleSnp2Arr[gpuId] >= _numSnp){
				_morePairsArr[gpuId] = false;
				break;
			}
		}

		return iter_block;
	}
	return 0;
}
