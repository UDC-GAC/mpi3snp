/*
 * Options.h
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_
#include "GPUMacros.h"

class Options {
public:
	Options();
	Options(int argc, char* argv[]);
	virtual ~Options();

	inline int getNumGPUs() {
		return _numGPUs;
	}

	inline uint16_t getGPUId(int id){

		if(id <= _numGPUs){
			return _GPUIds[id];
		}

		return _GPUIds[0];
	}

	inline string& getTPEDFileName() {
		return _tpedFileName;
	}

	inline string& getTFAMFileName(){
		return _tfamFileName;
	}

	inline string& getOutFileName() {
		return _outFileName;
	}

	inline uint16_t getNumOutputs(){
		return _numOutputs;
	}

	inline bool isHeteroGPUs(){
		return _heteroGPUs;
	}

	inline bool isMI(){
		return _filterMI;
	}

	// To parse the command line
	bool parse(int argc, char* argv[]);

	// To print the help for the users
	void printUsage();

private:
	string _tpedFileName;
	string _tfamFileName;
	string _outFileName;
	int _numGPUs;
	uint16_t* _GPUIds;
	uint16_t _numOutputs;
	bool _heteroGPUs;
	bool _filterMI;

	/*private member functions*/
	void _setDefaults();
};

#endif /* OPTIONS_H_ */
