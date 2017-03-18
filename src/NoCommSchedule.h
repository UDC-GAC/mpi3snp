//
// Created by christian on 16/03/17.
//

#ifndef MPI3SNP_NOCOMMSCHEDULE_H
#define MPI3SNP_NOCOMMSCHEDULE_H

#include "Block.h"
#include "Schedule.h"
#include <limits>

class NoCommSchedule {
public:
    NoCommSchedule(size_t problem_size, unsigned int k, unsigned int num_proc);

    ~NoCommSchedule();

    std::vector<std::vector<Block *>> Get_blocks(unsigned int proc_id);

    unsigned int Who_has(std::multiset<Block, std::less<Block>> c);

private:
    // Attributes
    static const unsigned int INVALID_PROC = std::numeric_limits<unsigned int>::max();
    const unsigned int k;
    std::vector<std::pair<std::multiset<Block, std::less<Block>>, unsigned int>> *dist;

    // Auxiliary methods
    void Create_blocks(size_t problem_size, unsigned int num_proc);

    void Combine_blocks(unsigned int k);

    void Generate_translation_table(int range, std::vector<std::pair<unsigned int, unsigned int>> *output);

    unsigned int Optimal_selection(std::vector<std::pair<unsigned int, unsigned int>> *tr, unsigned int range,
                                   unsigned int p, int flag);

    void Assign_combinations(unsigned int num_proc);
};


#endif //MPI3SNP_NOCOMMSCHEDULE_H
