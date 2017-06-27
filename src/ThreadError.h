//
// Created by christian on 27/06/17.
//

#ifndef MPI3SNP_THREADERROR_H
#define MPI3SNP_THREADERROR_H

class ThreadError : public std::runtime_error {
public:
    ThreadError(const std::string &message) :
            std::runtime_error(message) {};
};

#endif //MPI3SNP_THREADERROR_H
