/*
 * SearchMI.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 *
 *  Modified on: December 29, 2016
 *           by: Christian Ponte Fernández
 */

#ifndef SEARCHMI_H_
#define SEARCHMI_H_

#include "Search.h"
#include "ThreadParams.h"
#include "ContTable.h"
#include "float.h"
#include <mpi.h>
#include "IOMpi.h"

class SearchMI : public Search {
public:
    SearchMI(Options *options);

    virtual ~SearchMI();

    // Execute the epistasis search
    void execute();

private:
    static const int MPI_TAG_OUTPUT = 1;

    static void *_mpiMI(Options *options, vector<ThreadParams *> threadParams);

    static void *_threadMI(void *arg);

    vector<Engine *> _engines;
    vector<ThreadParams *> _threadParams;
};

#endif /* SEARCHMI_H_ */
