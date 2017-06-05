/*
 * GPUEngine.h
 *
 *  Created on: Oct 9, 2014
 *      Author: gonzales
 */

#ifndef GPUENGINE_H_
#define GPUENGINE_H_

#include "Options.h"
#include "GPUContTable.h"
#include "GPUVariables.h"
#include "MutualInfo.h"
#include "float.h"

class GPUEngine {
public:
	GPUEngine(Options *options);
	virtual ~GPUEngine();

	void setNums(uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls);

	// Initialize the engine
	void initialize(int gpuId);

	void loadSNPs(uint32_t* host0Cases, uint32_t* host1Cases, uint32_t* host2Cases,
			uint32_t* host0Ctrls, uint32_t* host1Ctrls, uint32_t* host2Ctrls);

	void mutualInfo(uint64_t numPairs, uint2* ids, MutualInfo* mutualInfo,
			float& minMI, uint16_t& minMIPos, uint16_t& numEntriesWithMI);

	inline uint16_t getGPUId(){
		return _gpuId;
	}

private:

	void _findNHighestMI(MutualInfo* mutualInfo, uint64_t numTotal,
		float& minMI, uint16_t& minMIPos,uint16_t& numEntriesWithMI);

	Options* _options;
	GPUInfo* _gpuInfo;
	uint16_t _gpuId;

	uint32_t _numSNPs;
	uint16_t _numCases;
	uint16_t _numCtrls;
	uint16_t _numEntriesCase;
	uint16_t _numEntriesCtrl;

	uint32_t* _dev0Cases;
	uint32_t* _dev1Cases;
	uint32_t* _dev2Cases;
	uint32_t* _dev0Ctrls;
	uint32_t* _dev1Ctrls;
	uint32_t* _dev2Ctrls;
	uint2* _devIds;

	// Auxiliary array for the contingency tables between the two kernels
	GPUDoubleContTable* _devDoubleTables;
	GPUDoubleContTable *_tables;

	// Auxiliary array to store the MI values of each block
	float *_devMIValues;
	float *_hostMIValues;

	// Auxiliary arrays to store the ids that are in the list of MIs
	uint3 *_devMiIds;
	uint3 *_hostMiIds;
};

#endif /* GPUENGINE_H_ */
