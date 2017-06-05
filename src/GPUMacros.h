/*
 * GPUMacros.h
 *
 *  Created on: Aug 31, 2012
 *      Author: yongchao
 */

#ifndef GPUMACROS_H_
#define GPUMACROS_H_
#include <cuda_runtime.h>
#include <cuda.h>
#include "Utils.h"

#define myCheckCudaError do{ \
	cudaError_t code = cudaGetLastError();\
	if(code != cudaSuccess) \
		Utils::log("Error %s occurred at line %d in function %s of file %s\n", cudaGetErrorString(code), __LINE__, __FUNCTION__, __FILE__);	\
	}while(0)

/*number of streaming processor scale*/
#define NUM_TH_PER_BLOCK	16

class GPUInfo
{
public:
	static GPUInfo* getGPUInfo() {
		if (_globalGPUs) {
			return _globalGPUs;
		}
		_globalGPUs = new GPUInfo();
		return _globalGPUs;

	}
	~GPUInfo() {
		_devices.clear();
	}
	inline int getNumGPUs() {
		return _devices.size();
	}
	inline int getGPUIndex(size_t i) {
		if (i >= _devices.size()) {
			Utils::exit(
					"GPU index (%d) is out of range (%d) in line %d in function %s\n",
					i, _devices.size(), __LINE__, __FUNCTION__);
		}
		return _devices[i].first;
	}

	/*select device*/
	inline void setDevice(int index) {
		if(index < 0 || index >= _devices.size()){
			Utils::exit("Incorrect GPU index (%d). Should be in [0, %d)\n", index, _devices.size());
		}

		index = _devices[index].first;
		cudaSetDevice(index);
		myCheckCudaError;
	}

	/*get GPU property*/
	inline cudaDeviceProp* getDeviceProp(int index)
	{
		return &(_devices[index].second);
	}

private:
	/*private member functions*/
	GPUInfo() {
		_devices.reserve(8);

		/*open the GPU device*/
		int numGPUs = 0;
		cudaGetDeviceCount(&numGPUs);
		myCheckCudaError;
		if (numGPUs == 0) {
			Utils::exit("There is no CUDA-enabled GPU available\n");
		} else {
			Utils::log("#CUDA-enabled devices: %d\n", numGPUs);
		}

		/*read the compatible device information*/
		cudaDeviceProp prop;
		for (int index = 0; index < numGPUs; ++index) {
			/*read the property of the device*/
			cudaGetDeviceProperties(&prop, index);
			myCheckCudaError;

			/*check if this GPU is compatible*/
#if defined(HAVE_SM_35)
			if ((prop.major * 10 + prop.minor) >= 35 && prop.canMapHostMemory){
#elif defined(HAVE_SM_30)
			if ((prop.major * 10 + prop.minor) >= 30 && prop.canMapHostMemory){
#else
			if (prop.major >= 2 && prop.canMapHostMemory){
#endif
				/*print out the GPU information*/
				Utils::log("---------Qualified GPU (index: %d)-------------\n", _devices.size());
				Utils::log("name:%s\n", prop.name);
				Utils::log("multiprocessor count:%d\n", prop.multiProcessorCount);
				Utils::log("clock rate:%d\n", prop.clockRate);
				Utils::log("shared memory:%d\n", prop.sharedMemPerBlock);
				Utils::log("global memory:%ld\n", prop.totalGlobalMem);
				Utils::log("registers per block:%d\n", prop.regsPerBlock);
				Utils::log("compute capability:%d.%d\n", prop.major, prop.minor);
				Utils::log("-----------------------------------------------\n");

				/*insert the qualified GPU*/
				insert(index, prop);
			}
		}
#if defined(HAVE_SM_35)
		Utils::log("Require GPUs with compute capability >= 3.5\n");
#elif defined(HAVE_SM_30)
		Utils::log("Require GPUs with compute capability >= 3.0\n");
#else
		Utils::log("Require GPU with compute capability >= 2.0\n");
#endif

		/*check the qualified GPUs*/
		if (_devices.size() == 0) {
			Utils::log(
					"NO qualified GPUs are available and will only use CPU\n");
		} else {
			Utils::log("Number of Qualified GPUs: %ld\n", _devices.size());
		}

		/*set the mapped memory support for each device*/
		for(size_t i = 0; i < _devices.size(); ++i){
			/*select device*/
			cudaSetDevice(_devices[i].first);
			myCheckCudaError;

		/*set the mapped memory support*/
			cudaSetDeviceFlags(cudaDeviceMapHost);
			myCheckCudaError;
		}
	}
	inline void insert(int gpuIndex, cudaDeviceProp& prop) {
		_devices.push_back(make_pair(gpuIndex, prop));
	}

	/*private member variables*/
	vector<pair<int, cudaDeviceProp> > _devices;
	static GPUInfo* _globalGPUs;
};

#endif /* GPUMACROS_H_ */
