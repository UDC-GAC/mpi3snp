/*
 * GPUSearchMI.h
 *
 *  Created on: Sep 30, 2014
 *      Author: gonzales
 */

#ifndef GPUSEARCHMI_H_
#define GPUSEARCHMI_H_

#include "Search.h"
#include <mpi.h>

class GPUSearchMI : public Search {
public:
    GPUSearchMI(Options *options);

    virtual ~GPUSearchMI();

    // Execute the epistasis search
    void execute();

private:
    MPI_Datatype MPI_MUTUAL_INFO;
};

#endif /* GPUSEARCHMI_H_ */
