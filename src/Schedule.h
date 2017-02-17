//
// Created by christian on 06/02/17.
//

#ifndef SCHEDULE_H
#define SCHEDULE_H

#include <vector>
#include "Block.h"

class Schedule {
public:
    virtual std::vector<std::vector<Block *>> getBlocks(int procId) =0;
    virtual int whoHas(std::vector<Block *> comb) =0;
};

#endif //SCHEDULE_H
