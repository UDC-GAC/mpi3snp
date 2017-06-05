/*
 * Options.cpp
 *
 *  Created on: 04/09/2014
 *      Author: jorge
 */

#include "Options.h"

Options::Options() {
	_setDefaults();
}

Options::Options(int argc, char* argv[]) {
	_setDefaults();

	// Parse the command line
	parse(argc, argv);
}

Options::~Options() {

	if(_GPUIds){
		delete [] _GPUIds;
	}

}

void Options::printUsage() {
	fprintf(stderr,
			"--------------------------------------------------------------------------------\n\n");
	fprintf(stderr,
			"GPU3SNP (v%s) is a tool to detect 3-way interaction in Genome Wide Association datasets using"
			"GPUs\n",	VERSION);

	fprintf(stderr,
			"--------------------------------------------------------------------------------\n\n");

	fprintf(stderr, "Usage: GPU3SNP -rp srcTPEDFile -rf srcTFAMFile -o outFile [options]\n");

	// Input
	fprintf(stderr, "\tsrcTPEDFile: the file name for the SNPs inputs with -tped format\n");
	fprintf(stderr, "\tsrcTFAMFile: the file name for the individuals information with -tfam format\n");

	// Output
	fprintf(stderr, "Output:\n");
	fprintf(stderr, "\toutFile: output file name\n");

	// Compute
	fprintf(stderr, "Computational options:\n");
	fprintf(stderr, "\t-t <int> (number of GPUs for computation, default = %d)\n",
			_numGPUs);

	fprintf(stderr, "\t-g <int> [int] ... (index of the GPUs. There must be the same number of values as GPUs indicated with -t. Default = %d)\n",
			_GPUIds[0]);

	fprintf(stderr, "\t-no <int> (number of outputs, default = %hu)\n",
			_numOutputs);

	fprintf(stderr, "\t-het (to indicate that there are GPUs of different types and a dynamic distribution is applied)\n");

	fprintf(stderr, "\t-ig (apply Information Gain for filtering)\n");

	fprintf(stderr, "\t-mi (apply Mutual Information for filtering, by default)\n");

	fprintf(stderr, "Others:\n");
	fprintf(stderr, "\t-h <print out the usage of the program)\n");
}

void Options::_setDefaults() {
	// Empty string means outputing to STDOUT
	_tpedFileName = "";
	_tfamFileName = "";
	_outFileName = "";
	_numGPUs = 1;

	_GPUIds = new uint16_t[1];
	_GPUIds[0] = 0;

	_numOutputs = 10;
	_heteroGPUs = false;
	_filterMI = true;
}

bool Options::parse(int argc, char* argv[]) {
	int intVal, i, j;
	uint16_t uintVal;
	int argind = 1;

	if (argc < 2) {
		return false;
	}

	// Check the help
	if (!strcmp(argv[argind], "-h") || !strcmp(argv[argind], "-?")) {
		return false;
	}

	// Print out the command line
	fprintf(stderr, "Command: ");
	for (int i = 0; i < argc; ++i) {
		fprintf(stderr, "%s ", argv[i]);
	}
	fputc('\n', stderr);

	/*check the availability of GPUs*/
	GPUInfo *gpuInfo = GPUInfo::getGPUInfo();
	if (gpuInfo->getNumGPUs() == 0) {
		Utils::exit("No compatible GPUs are available in your machine\n");
	}

	// For the other options
	while (argind < argc) {
		// Input file
		if (!strcmp(argv[argind], "-rp")) {
			argind++;
			if (argind < argc) {
				_tpedFileName = argv[argind];
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
				return false;
			}
		}
		else if (!strcmp(argv[argind], "-rf")) {
			argind++;
			if (argind < argc) {
				_tfamFileName = argv[argind];
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
				return false;
			}
		}
		else if (!strcmp(argv[argind], "-o")) {
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_outFileName = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-t")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_numGPUs = intVal;
				delete [] _GPUIds;
				_GPUIds = new uint16_t[_numGPUs];
			} else {
				Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-g")) {
			argind++;
			for(i=0; i<_numGPUs; i++){
				if (argind < argc) {
					sscanf(argv[argind], "%d", &intVal);
					if ((intVal < 0) || (intVal>=gpuInfo->getNumGPUs())){
						Utils::log("value %d not valid for GPU index\n", intVal);
						return false;
					}
					for(j=0; j<i; j++){
						if(_GPUIds[j] == intVal){
							Utils::log("repeating GPU indexes not allowed\n");
							return false;
						}
					}
					_GPUIds[i] = intVal;
				} else {
					Utils::log("not specify enough values for the parameter -g\n");
					return false;
				}
				argind++;
			}
		}
		else if (!strcmp(argv[argind], "-no")) {
			uintVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%hu", &uintVal);
				argind++;
				_numOutputs = uintVal;
			} else {
				Utils::log("not specify value for the parameter %s\n", argv[argind - 1]);
				return false;
			}
		}
		else if (!strcmp(argv[argind], "-het")) {
			_heteroGPUs = true;
			argind++;
		}
		else if (!strcmp(argv[argind], "-ig")) {
			_filterMI = false;
			argind++;
		}
		else if (!strcmp(argv[argind], "-mi")) {
			_filterMI = true;
			argind++;
		}
	}

	if(_filterMI){
		Utils::log("Applying Mutual Information\n");
	}
	else{
		Utils::log("Applying Information Gain\n");
	}
	Utils::log("Number of GPUs: %d\n", _numGPUs);
	Utils::log("GPU index:");
	for(i=0; i<_numGPUs; i++){
		Utils::log(" %d", _GPUIds[i]);
	}
	Utils::log("\nNumber of outputs: %hu\n", _numOutputs);
	Utils::log("Number of pairs by block: %lu\n", NUM_PAIRS_BLOCK);
	if(_numGPUs > 1){
		if(_heteroGPUs){
			Utils::log("Dynamic distribution among GPUs\n");
		}
		else{
			Utils::log("Static distribution among GPUs\n");
		}
	}

	if(!strcmp(_tfamFileName.c_str(), "") || !strcmp(_tpedFileName.c_str(), "")){
		Utils::log("Input files not specified!!!\n");
		return false;
	}

	if(!strcmp(_outFileName.c_str(), "")){
		Utils::log("Output file not specified!!!\n");
		return false;
	}

	return true;
}
