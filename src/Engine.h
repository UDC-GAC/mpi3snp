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
    class Error : public std::runtime_error {
    public:
        explicit Error(const std::string &message) : runtime_error(message) {};

        ~Error() override = default;
    };

    virtual void run(std::string tped, std::string tfam, std::vector<MutualInfo> &mutual_info, size_t num_outputs,
                     Statistics &statistics) = 0;
};

#endif //MPI3SNP_ENGINE_H
