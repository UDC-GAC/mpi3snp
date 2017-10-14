//
// Created by christian on 13/10/17.
//

#ifndef MPI3SNP_ENGINE_H
#define MPI3SNP_ENGINE_H


#include <vector>
#include "MutualInfo.h"
#include "Statistics.h"

class Engine {
public:
    virtual void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs,
                     Statistics &statistics) = 0;
};

#endif //MPI3SNP_ENGINE_H
