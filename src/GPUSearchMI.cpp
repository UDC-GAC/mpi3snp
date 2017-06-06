/*
 * GPUSearchMI.cpp
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#include "GPUSearchMI.h"
#include "GPUEngine.h"

GPUSearchMI::GPUSearchMI(Options *options) : Search(options) {}

GPUSearchMI::~GPUSearchMI() {}

void GPUSearchMI::execute() {
    GPUEngine gpu_engine(_options);
    gpu_engine.run();
}