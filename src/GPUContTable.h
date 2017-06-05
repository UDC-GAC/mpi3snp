/*
 * GPUContTable.h
 *
 *  Created on: 09/10/2014
 *      Author: jorge
 */

#ifndef CONTTABLE_H_
#define CONTTABLE_H_

#include "Macros.h"

/*
 * Structure for the auxiliar 2-SNP contingency tables
 */

struct GPUDoubleContTable
{
	GPUDoubleContTable(){

	}

	~GPUDoubleContTable(){
	}

	void initialize(uint16_t numEntriesCase, uint16_t numEntriesCtrl){

		cudaMalloc(&_cases00, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases01, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases02, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases10, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases11, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases12, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases20, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases21, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_cases22, numEntriesCase*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls00, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls01, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls02, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls10, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls11, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls12, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls20, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls21, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;

		cudaMalloc(&_ctrls22, numEntriesCtrl*sizeof(uint32_t));
		myCheckCudaError;
	}

	void finalize(){
		if(_cases00){
			cudaFree(_cases00);
			myCheckCudaError;
		}
		if(_cases01){
			cudaFree(_cases01);
			myCheckCudaError;
		}
		if(_cases02){
			cudaFree(_cases02);
			myCheckCudaError;
		}
		if(_cases10){
			cudaFree(_cases10);
			myCheckCudaError;
		}
		if(_cases11){
			cudaFree(_cases11);
			myCheckCudaError;
		}
		if(_cases12){
			cudaFree(_cases12);
			myCheckCudaError;
		}
		if(_cases20){
			cudaFree(_cases20);
			myCheckCudaError;
		}
		if(_cases21){
			cudaFree(_cases21);
			myCheckCudaError;
		}
		if(_cases22){
			cudaFree(_cases22);
			myCheckCudaError;
		}
		if(_ctrls00){
			cudaFree(_ctrls00);
			myCheckCudaError;
		}
		if(_ctrls01){
			cudaFree(_ctrls01);
			myCheckCudaError;
		}
		if(_ctrls02){
			cudaFree(_ctrls02);
			myCheckCudaError;
		}
		if(_ctrls10){
			cudaFree(_ctrls10);
			myCheckCudaError;
		}
		if(_ctrls11){
			cudaFree(_ctrls11);
			myCheckCudaError;
		}
		if(_ctrls12){
			cudaFree(_ctrls12);
			myCheckCudaError;
		}
		if(_ctrls20){
			cudaFree(_ctrls20);
			myCheckCudaError;
		}
		if(_ctrls21){
			cudaFree(_ctrls21);
			myCheckCudaError;
		}
		if(_ctrls22){
			cudaFree(_ctrls22);
			myCheckCudaError;
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


#endif /* CONTTABLE_H_ */
