/*
 * SearchMI.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#include "SearchMI.h"

SearchMI::SearchMI(Options* options) : Search(options){
	_engines.resize(options->getNumCPUs());
	_threadParams.resize(options->getNumCPUs());

	for (int tid = 0; tid < options->getNumCPUs(); ++tid) {
		_engines[tid] = new Engine(options);
		// Create parameters for CPU threads
		_threadParams[tid] = new ThreadParams(tid, _engines[tid]);
	}
}

SearchMI::~SearchMI() {
	for(int i=0; i<_engines.size(); i++){
		delete _engines[i];
		delete _threadParams[i];
	}

	_engines.clear();
	_threadParams.clear();
}

void SearchMI::execute(){
	double stime = Utils::getSysTime();
	double etime;

	SNPDistributor* distributor = new SNPDistributor(_options);
	distributor->loadSNPSet();

	etime = Utils::getSysTime();
	Utils::log("Loaded %ld SNPs in %.2f seconds\n", distributor->getNumSnp(), etime - stime);

	vector<pthread_t> threadIDs(_options->getNumCPUs(), 0);

	// Computation of the single-SNP entropy
	for(int tid=0; tid<_options->getNumCPUs(); tid++) {
		// All threads use the same distributor so it gives them the correct pair of SNPs
		_threadParams[tid]->init(distributor, _options->getNumOutputs());

		// Create thread entities that call to the functions below
		if (pthread_create(&threadIDs[tid], NULL, _threadMI, _threadParams[tid]) != 0) {
			Utils::exit("Thread creating failed\n");
		}
	}

	MutualInfo* auxMutualInfo = new MutualInfo[_options->getNumCPUs()*_options->getNumOutputs()];

	// Wait for the completion of all threads
	for(int tid=0; tid<_options->getNumCPUs(); tid++) {
		pthread_join(threadIDs[tid], NULL);
		memcpy(&auxMutualInfo[tid*_options->getNumOutputs()], _threadParams[tid]->_mutualInfo,
				_options->getNumOutputs()*sizeof(MutualInfo));
	}

	// Sort the auxiliar array and print the results
	std::sort(auxMutualInfo, auxMutualInfo+_options->getNumOutputs()*_options->getNumCPUs());
	distributor->printMI(auxMutualInfo+_options->getNumOutputs()*(_options->getNumCPUs()-1),
			_options->getNumOutputs());

	Utils::log("3-SNP analysis finalized\n");

#ifdef DEBUG
	uint32_t numAnalyzed = 0;

	for(int tid=0; tid<_options->getNumCPUs(); tid++){
		numAnalyzed+=_threadParams[tid]->_numAnalyzed;
	}

	Utils::log("Total analysis: %" PRIu64 "\n", numAnalyzed);
#endif

	// Release the distributor
	delete distributor;
	delete [] auxMutualInfo;
	threadIDs.clear();
}

void* SearchMI::_threadMI(void* arg){
	ThreadParams* params = (ThreadParams*) arg;

	Engine* engine = (Engine*) params->_engine;

	SNPDistributor* distributor = params->_distributor;
	uint32_t numSnp = distributor->getNumSnp();
	uint16_t numCases = distributor->getNumCases();
	uint16_t numCtrls = distributor->getNumCtrls();
	uint16_t numOutputs = params->_numOutputs;
	engine->setNums(numSnp, numCases, numCtrls);

	vector<SNP*> snpSet = distributor->getSnpSet();
	// In this case the ids are necessary to access to the single-SNP entropy
	uint32_t* auxIds = new uint32_t[2*NUM_PAIRS_BLOCK];

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
		numPairsBlock = distributor->getPairsSNPs(auxIds);

		if(numPairsBlock){
			myTotalAnal += engine->mutualInfo(snpSet, numPairsBlock, auxIds, mutualInfo, numOutputs, minMI,
					minMIPos, numEntriesWithMI);
		}
		else{
			moreAnal = false;
		}
	}

	params->_numAnalyzed = myTotalAnal;
	memcpy(params->_mutualInfo, mutualInfo, numOutputs*sizeof(MutualInfo));

#ifdef BENCHMARKING
	etime = Utils::getSysTime();
	Utils::log("CPU thread (%d) %f seconds calculating %lu analysis\n",
			params->_tid, etime-stime, myTotalAnal);
#endif

	delete [] auxIds;
	delete [] mutualInfo;
	return NULL;
}
