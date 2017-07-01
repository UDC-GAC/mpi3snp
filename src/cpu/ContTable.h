/*
 * ContTable.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef CONTTABLE_H_
#define CONTTABLE_H_

#include <cstring>

/*
 * Structure for the auxiliar 2-SNP contingency tables
 */

struct DoubleContTable
{
	DoubleContTable(uint16_t numCases, uint16_t numCtrls){
		uint16_t numEntryCases = numCases/32+((numCases%32)>0);
		uint16_t numEntryCtrls = numCtrls/32+((numCtrls%32)>0);

		_cases00 = new uint32_t[numEntryCases];
		if (_cases00 == NULL) {
			exit(1);
		}

		_cases01 = new uint32_t[numEntryCases];
		if (_cases01 == NULL) {
			exit(1);
		}

		_cases02 = new uint32_t[numEntryCases];
		if (_cases02 == NULL) {
			exit(1);
		}

		_cases10 = new uint32_t[numEntryCases];
		if (_cases10 == NULL) {
			exit(1);
		}

		_cases11 = new uint32_t[numEntryCases];
		if (_cases11 == NULL) {
			exit(1);
		}

		_cases12 = new uint32_t[numEntryCases];
		if (_cases12 == NULL) {
			exit(1);
		}

		_cases20 = new uint32_t[numEntryCases];
		if (_cases20 == NULL) {
			exit(1);
		}

		_cases21 = new uint32_t[numEntryCases];
		if (_cases21 == NULL) {
			exit(1);
		}

		_cases22 = new uint32_t[numEntryCases];
		if (_cases22 == NULL) {
			exit(1);
		}

		_ctrls00 = new uint32_t[numEntryCtrls];
		if (_ctrls00 == NULL) {
			exit(1);
		}

		_ctrls01 = new uint32_t[numEntryCtrls];
		if (_ctrls01 == NULL) {
			exit(1);
		}

		_ctrls02 = new uint32_t[numEntryCtrls];
		if (_ctrls02 == NULL) {
			exit(1);
		}

		_ctrls10 = new uint32_t[numEntryCtrls];
		if (_ctrls10 == NULL) {
			exit(1);
		}

		_ctrls11 = new uint32_t[numEntryCtrls];
		if (_ctrls11 == NULL) {
			exit(1);
		}

		_ctrls12 = new uint32_t[numEntryCtrls];
		if (_ctrls12 == NULL) {
			exit(1);
		}

		_ctrls20 = new uint32_t[numEntryCtrls];
		if (_ctrls20 == NULL) {
			exit(1);
		}

		_ctrls21 = new uint32_t[numEntryCtrls];
		if (_ctrls21 == NULL) {
			exit(1);
		}

		_ctrls22 = new uint32_t[numEntryCtrls];
		if (_ctrls22 == NULL) {
			exit(1);
		}
	}

	uint32_t* _cases00;
	uint32_t* _cases01;
	uint32_t* _cases02;
	uint32_t* _cases10;
	uint32_t* _cases11;
	uint32_t* _cases12;
	uint32_t* _cases20;
	uint32_t* _cases21;
	uint32_t* _cases22;
	uint32_t* _ctrls00;
	uint32_t* _ctrls01;
	uint32_t* _ctrls02;
	uint32_t* _ctrls10;
	uint32_t* _ctrls11;
	uint32_t* _ctrls12;
	uint32_t* _ctrls20;
	uint32_t* _ctrls21;
	uint32_t* _ctrls22;
};

struct TripleContTable
{
	TripleContTable(){
		clearValues();
	}

	void clearValues(){
		int i;
		for(i=0; i<27; i++){
			_cases[i] = 0;
		}
		for(i=0; i<27; i++){
			_ctrls[i] = 0;
		}
	}

	void setValues(uint16_t* cases, uint16_t* ctrls){
		std::memcpy(_cases, cases, 27*sizeof(uint16_t));
		std::memcpy(_ctrls, ctrls, 27*sizeof(uint16_t));
	}

	uint16_t _cases[27];
	uint16_t _ctrls[27];
};


#endif /* CONTTABLE_H_ */
