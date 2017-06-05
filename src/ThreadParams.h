/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "Macros.h"
#include "GPUSNPDistributorStatic.h"
#include "GPUEngine.h"
#include "MutualInfo.h"

struct ThreadParams
{
	ThreadParams(int tid, GPUEngine* engine, bool dynamicDist){
		_tid = tid;
		_distributor = NULL;
		_engine = engine;
		_numOutputs = 0;
		_numAnalyzed = 0;
		_mutualInfo = NULL;
		_dynamicDist = dynamicDist;
	}

	~ThreadParams()
	{
		delete [] _mutualInfo;
	}

	inline void init(GPUSNPDistributor* distributor, uint16_t numOutputs)
	{
		_distributor = distributor;
		_numOutputs = numOutputs;
		_mutualInfo = new MutualInfo[numOutputs];
	}

	int _tid;
	GPUSNPDistributor* _distributor;
	GPUEngine* _engine;

	// Values for the return information
	uint16_t _numOutputs;
	uint64_t _numAnalyzed;
	MutualInfo* _mutualInfo;
	bool _dynamicDist;
};


#endif /* THREADPARAMS_H_ */
