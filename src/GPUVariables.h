/*
 * GPUVariables.h
 *
 *  Created on: Oct 09, 2014
 *      Author: gonzales
 */

#ifndef GPUVARIABLES_H_
#define GPUVARIABLES_H_

// Constant variable to keep the information of the entropy of the cases/controls
extern __constant__ float _invInds;
extern __constant__ float _entY;
extern __constant__ float _MAX_FLOAT;

#endif /* GPUVARIABLES_H_ */
