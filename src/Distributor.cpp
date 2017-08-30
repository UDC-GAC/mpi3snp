//
// Created by christian on 16/06/17.
//

#include "Distributor.h"
#include <cstdlib>

Distributor::Distributor(unsigned int p_num, unsigned int p_id, unsigned long count) :
        count(count) {
    unsigned long dist_size = count / p_num + (p_id < (count % p_num));
    distribution.resize(dist_size);
    distribution[0] = p_id;
    distribution[1] = 2 * p_num - p_id - 1;
    for (unsigned int i = 2; i < dist_size; i++) {
        distribution[i] = distribution[i - 2] + 2 * p_num;
    }
}

void Distributor::Get_pairs(unsigned int t_num, unsigned int t_id,
                                     std::vector<std::pair<uint32_t, uint32_t>> &values) {
    for (size_t i = t_id; i < distribution.size(); i += t_num) {
        for (uint32_t j = distribution[i] + 1; j < count - 1; j++) {
            values.push_back(std::make_pair(distribution[i], j));
        }
    }
}