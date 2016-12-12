/*
 * ThreadParam.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef THREADPARAMS_H_
#define THREADPARAMS_H_

#include "Macros.h"
#include "SNPDistributor.h"
#include "Engine.h"
#include "MutualInfo.h"

struct ThreadParams
{
	ThreadParams(int tid, Engine* engine){
		_tid = tid;
		_distributor = NULL;
		_engine = engine;
		_numOutputs = 0;
		_numAnalyzed = 0;
		_mutualInfo = NULL;
	}

	~ThreadParams()
	{
		delete _mutualInfo;
	}

	inline void init(SNPDistributor* distributor, uint16_t numOutputs)
	{
		_distributor = distributor;
		_numOutputs = numOutputs;
		_mutualInfo = new MutualInfo[numOutputs];
	}

	int _tid;
	SNPDistributor* _distributor;
	Engine* _engine;

	// Values for the return information
	uint16_t _numOutputs;
	uint64_t _numAnalyzed;
	MutualInfo* _mutualInfo;
};


#endif /* THREADPARAMS_H_ */
