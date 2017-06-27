/*
 * Engine.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include "Macros.h"
#include "../SNP.h"
#include "Options.h"
#include "../MutualInfo.h"
#include "ContTable.h"

class Engine {
public:
	Engine(Options* options);
	virtual ~Engine();

	inline void setNums(uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls){
		_numSNPs = numSNPs;
		_numCases = numCases;
		_numCtrls = numCtrls;
		_numEntriesCases = numCases/32+((numCases%32)>0);
		_numEntriesCtrls = numCtrls/32+((numCtrls%32)>0);
		_invInds = 1.0/(numCases+numCtrls);

		float p = numCases*_invInds;
		_entY = (-1.0)*p*log2(p);

		p = numCtrls*_invInds;
		_entY -= p*log2(p);
	}

	// This function works with the whole SNPSet in memory
	uint64_t mutualInfo(vector<SNP*> snpSet, uint64_t numPairs, uint32_t* ids, MutualInfo* mutualInfo,
			uint16_t numOutputs, float& minMI, uint16_t& minMIPos, uint16_t& numEntriesWithMI);

private:

	void _fillDoubleContTable(SNP* snp1, SNP* snp2, DoubleContTable* table);

	float _calcDoubleMI(DoubleContTable* table);

	void _fillTripleContTable(DoubleContTable* doubleTable, TripleContTable* tripleTable, SNP* snp3);

	float _calcTripleMI(TripleContTable* table);

	Options* _options;

	uint32_t _numSNPs;
	uint16_t _numCases;
	uint16_t _numCtrls;
	uint16_t _numEntriesCases;
	uint16_t _numEntriesCtrls;
	float _invInds;

	// Entropy of Y
	float _entY;
};

#endif /* ENGINE_H_ */
