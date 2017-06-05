/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include "GPUSearchMI.h"

GPUSearchMI::GPUSearchMI(Options* options) : Search(options){
	_engines.resize(options->getNumGPUs());
	_threadParams.resize(options->getNumGPUs());
	_dynamicDist = options->isHeteroGPUs();

	for (int tid = 0; tid < options->getNumGPUs(); ++tid) {
		_engines[tid] = new GPUEngine(options);
		// Create parameters for CPU threads
		_threadParams[tid] = new ThreadParams(tid, _engines[tid], _dynamicDist);
	}
}

GPUSearchMI::~GPUSearchMI() {
	for(int i=0; i<_engines.size(); i++){
		delete _engines[i];
		delete _threadParams[i];
	}

	_engines.clear();
	_threadParams.clear();
}

void GPUSearchMI::execute(){
	double stime = Utils::getSysTime();
	double etime;

	GPUSNPDistributor* distributor;

	if(_dynamicDist){
		distributor = new GPUSNPDistributor(_options);
	}
	else{
		distributor = new GPUSNPDistributorStatic(_options);
	}

	distributor->loadSNPSet();

	etime = Utils::getSysTime();
	Utils::log("Loaded %ld SNPs in %.2f seconds\n", distributor->getNumSnp(), etime - stime);

	vector<pthread_t> threadIDs(_options->getNumGPUs(), 0);

	for(int tid=0; tid<_options->getNumGPUs(); tid++) {
		// All threads use the same distributor so it gives them the correct pair of SNPs
		_threadParams[tid]->init(distributor, _options->getNumOutputs());

		// Create thread entities that call to the functions below
		if (pthread_create(&threadIDs[tid], NULL, _GPUMI, _threadParams[tid]) != 0) {
			Utils::exit("Thread creating failed\n");
		}
	}

	MutualInfo* auxMutualInfo = new MutualInfo[_options->getNumGPUs()*_options->getNumOutputs()];

	// Wait for the completion of all threads
	for(int tid=0; tid<_options->getNumGPUs(); tid++) {
		pthread_join(threadIDs[tid], NULL);
		memcpy(&auxMutualInfo[tid*_options->getNumOutputs()], _threadParams[tid]->_mutualInfo,
				_options->getNumOutputs()*sizeof(MutualInfo));
	}

	// Sort the auxiliar array and print the results
	std::sort(auxMutualInfo, auxMutualInfo+_options->getNumOutputs()*_options->getNumGPUs());
	distributor->printMI(auxMutualInfo+_options->getNumOutputs()*(_options->getNumGPUs()-1),
			_options->getNumOutputs());

	Utils::log("3-SNP analysis finalized\n");

#ifdef DEBUG
	uint32_t numAnalyzed = 0;

	for(int tid=0; tid<_options->getNumGPUs(); tid++){
		numAnalyzed+=_threadParams[tid]->_numAnalyzed;
	}

	Utils::log("Total analysis: %" PRIu64 "\n", numAnalyzed);
#endif

	// Release the distributor
	delete distributor;
	delete [] auxMutualInfo;
	threadIDs.clear();
}

void* GPUSearchMI::_GPUMI(void* arg){
	ThreadParams* params = (ThreadParams*) arg;

	GPUSNPDistributor* distributor;
	if(params->_dynamicDist){
		distributor = (GPUSNPDistributor *)params->_distributor;
	}
	else{
		distributor = (GPUSNPDistributorStatic *)params->_distributor;
	}
	GPUEngine* engine = (GPUEngine*) params->_engine;
	//initialize GPU
	engine->initialize(params->_tid);
	uint32_t numSnp = distributor->getNumSnp();
	uint16_t numCases = distributor->getNumCases();
	uint16_t numCtrls = distributor->getNumCtrls();
	uint16_t numOutputs = params->_numOutputs;
	engine->setNums(numSnp, numCases, numCtrls);
	engine->loadSNPs(distributor->getHost0Cases(), distributor->getHost1Cases(),
			distributor->getHost2Cases(), distributor->getHost0Ctrls(),
			distributor->getHost1Ctrls(), distributor->getHost2Ctrls());

	uint2* auxIds;
	cudaMallocHost(&auxIds, NUM_PAIRS_BLOCK*sizeof(uint2));
	myCheckCudaError;

	// Variables to work with the outputs
	MutualInfo* mutualInfo = new MutualInfo[numOutputs];
	// The minimum value in the array
	float minMI = FLT_MAX;
	// The position of the minimum value
	uint16_t minMIPos = 0;
	// Number of entries of the array full
	uint16_t numEntriesWithMI = 0;

	bool moreAnal = true;
	uint64_t myTotalAnal = 0;
	uint64_t numPairsBlock = 0;

#ifdef BENCHMARKING
	double stime = Utils::getSysTime();
	double etime;
#endif

	while(moreAnal){
		// Take some SNPs
		numPairsBlock = distributor->getPairsSNPs(auxIds, myTotalAnal, params->_tid);

		if(numPairsBlock){
			engine->mutualInfo(numPairsBlock, auxIds, mutualInfo, minMI,
					minMIPos, numEntriesWithMI);
		}
		else{
			moreAnal = false;
		}
	}

#ifdef BENCHMARKING
	cudaDeviceSynchronize();
	etime = Utils::getSysTime();
	Utils::log("GPU thread (%d) %f seconds calculating %lu analysis\n",
			params->_tid, etime-stime, myTotalAnal);
#endif

	params->_numAnalyzed = myTotalAnal;
	memcpy(params->_mutualInfo, mutualInfo, numOutputs*sizeof(MutualInfo));

	cudaFreeHost(auxIds);
	myCheckCudaError;
	delete [] mutualInfo;
	return NULL;
}
