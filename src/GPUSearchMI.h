/*
 * GPUSearchMI.h
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#ifndef GPUSEARCHMI_H_
#define GPUSEARCHMI_H_

#include "Search.h"

class GPUSearchMI : public Search {
public:
	GPUSearchMI(Options* options);
	virtual ~GPUSearchMI();

	// Execute the epistasis search
	void execute();
};

#endif /* GPUSEARCHMI_H_ */
