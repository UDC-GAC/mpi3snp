//
// Created by christian on 16/06/17.
//

#include "Distributor.h"

Distributor::Distributor(unsigned int proc_num, unsigned int proc_id, unsigned long item_num, bool shared,
                         unsigned long block_size) :
        item_num(item_num),
        items_left(true),
        block_size(block_size),
        shared(shared) {
    pthread_mutex_init(&lock, nullptr);

    dist_size = item_num / proc_num +
                (proc_id < (item_num % proc_num));
    dist = new uint32_t[dist_size];
    dist[0] = proc_id;
    dist[1] = 2 * proc_num - proc_id - 1;
    for (int i = 2; i < dist_size; i++) {
        dist[i] = dist[i - 2] + 2 * proc_num;
    }

    dist_it = 0;
    it1 = dist[dist_it++];
    it2 = it1 + 1;
}

Distributor::~Distributor() {
    delete[] dist;
    pthread_mutex_destroy(&lock);
}

unsigned long Distributor::Get_pairs(std::vector<std::pair<uint32_t,uint32_t>> &values, uint64_t &totalAnal) {
    if (shared) {
        return Get_pairs_lock(values, totalAnal);
    } else {
        return Get_pairs_nolock(values, totalAnal);
    }
}

unsigned long Distributor::Get_pairs_nolock(std::vector<std::pair<uint32_t,uint32_t>> &values, uint64_t &totalAnal) {
    unsigned long iter_block = 0;

    if (items_left) {
        while (dist_it < dist_size) {
            while (it2 < item_num - 1 && iter_block < block_size) {
                totalAnal += item_num - it2 - 1;
                values.push_back(std::make_pair(it1, it2++));
                iter_block++;
            }

            if (iter_block == block_size) {
                return iter_block;
            } else {
                it1 = dist[dist_it++];
                it2 = it1 + 1;
            }
        }
        items_left = false;
    }

    return iter_block;
}

unsigned long Distributor::Get_pairs_lock(std::vector<std::pair<uint32_t,uint32_t>> &values, uint64_t &totalAnal) {
    unsigned long iter_block = 0;

    pthread_mutex_lock(&lock);
    iter_block = Get_pairs_nolock(values, totalAnal);
    pthread_mutex_unlock(&lock);

    return iter_block;
}