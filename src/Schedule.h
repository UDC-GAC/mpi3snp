//
// Created by christian on 06/02/17.
//

#ifndef SCHEDULE_H
#define SCHEDULE_H

#include <vector>
#include <set>
#include "Block.h"

class Schedule {
public:
    virtual void Get_blocks(unsigned int proc_id, std::vector<std::multiset<Block, std::less<Block>>> *block_list) =0;
    virtual unsigned int Who_has(std::multiset<Block, std::less<Block>> c) =0;
    virtual ~Schedule(){};
};

#endif //SCHEDULE_H
