//
// Created by christian on 16/06/17.
//

#ifndef MPI3SNP_DISTRIBUTOR_H
#define MPI3SNP_DISTRIBUTOR_H

#include <inttypes.h>
#include <pthread.h>
#include <vector>

class Distributor {
public:
    Distributor(unsigned int proc_num, unsigned int proc_id, unsigned long item_num, bool shared,
                unsigned long block_size = DEFAULT_PAIRS_BLOCK);

    ~Distributor();

    unsigned long Get_pairs(std::vector<std::pair<uint32_t,uint32_t>> &values, uint64_t &totalAnal);

    static const unsigned long DEFAULT_PAIRS_BLOCK = 25000;

private:
    unsigned long Get_pairs_lock(std::vector<std::pair<uint32_t,uint32_t>> &values, uint64_t &totalAnal);

    unsigned long Get_pairs_nolock(std::vector<std::pair<uint32_t,uint32_t>> &values, uint64_t &totalAnal);

    bool items_left;
    unsigned long item_num;
    unsigned long block_size;
    bool shared;
    pthread_mutex_t lock;

    // Iterators for the SNPs
    unsigned long dist_size, dist_it;
    uint32_t *dist;
    uint32_t it1, it2;
};


#endif //MPI3SNP_DISTRIBUTOR_H
