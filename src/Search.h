/*
 * Search.h
 *
 *  Created on: Sep 4, 2014
 *      Author: gonzales
 */

#ifndef SEARCH_H_
#define SEARCH_H_

#include "Options.h"

class Search {
public:
	Search(Options* options);
	virtual ~Search();

	// Execute the epistasis search
	virtual void execute() = 0;

protected:
	Options* _options;
};

#endif /* SEARCH_H_ */
