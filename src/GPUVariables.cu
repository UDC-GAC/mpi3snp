/*
 * GPUVariables.cu
 *
 *  Created on: Oct 09, 2014
 *      Author: gonzales
 */

#include "GPUMacros.h"
#include "GPUVariables.h"

// Constant variable to keep the information of the entropy of the cases/controls
__constant__ float _invInds;
__constant__ float _entY;
__constant__ float _MAX_FLOAT;

