/*
 * SearchMI.h
 *
 *  Created on: Sep 5, 2014
 *      Author: gonzales
 */

#ifndef SEARCHMI_H_
#define SEARCHMI_H_

#include "Search.h"
#include "ThreadParams.h"
#include "ContTable.h"
#include "float.h"

class SearchMI : public Search {
public:
	SearchMI(Options* options);
	virtual ~SearchMI();

	// Execute the epistasis search
	void execute();

private:
	static void* _threadMI(void* arg);

	vector<Engine*> _engines;
	vector<ThreadParams*> _threadParams;
};

#endif /* SEARCHMI_H_ */
