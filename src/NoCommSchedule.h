//
// Created by christian on 16/03/17.
//

#ifndef MPI3SNP_NOCOMMSCHEDULE_H
#define MPI3SNP_NOCOMMSCHEDULE_H

#include "Block.h"
#include "Schedule.h"

class NoCommSchedule {
public:
    NoCommSchedule(size_t problem_size, unsigned int k, unsigned int num_proc);

    ~NoCommSchedule();

    std::vector<std::vector<Block *>> Get_blocks(unsigned int proc_id);

    int Who_has(std::vector<Block *> comb);

    unsigned int Who_has2(std::multiset<Block, std::less<Block>> c);

private:
    unsigned int k;

    // Store block distribution
    std::vector<std::pair<std::multiset<Block, std::less<Block>>, unsigned int>> *dist;

    // Auxiliary methods
    void distribute_problem(size_t problem_size, unsigned int num_proc);

    void combine_blocks(unsigned int k);

    void assign_combinations();

    bool Multiset_greater_than(std::multiset<Block> m1, std::multiset<Block> m2);
};


#endif //MPI3SNP_NOCOMMSCHEDULE_H
