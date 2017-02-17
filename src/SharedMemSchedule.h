//
// Created by christian on 14/02/17.
//

#ifndef SHAREDMEMSCHEDULE_H
#define SHAREDMEMSCHEDULE_H

#include "Schedule.h"
#include <cstdint>

class SharedMemSchedule : public Schedule {
public:
    SharedMemSchedule(size_t size, int numProc);

    ~SharedMemSchedule();

    std::vector<std::vector<Block *>> getBlocks(int procId);

    int whoHas(std::vector<Block *> comb);

private:
    size_t size;
    int numProc;

    // Store block distribution
    std::vector<Block *> blockDiv;
    std::vector<std::vector<Block *>> comb;
    std::vector<int> proc;
};


#endif //SHAREDMEMSCHEDULE_H
