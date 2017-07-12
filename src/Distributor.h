//
// Created by christian on 16/06/17.
//

#ifndef MPI3SNP_DISTRIBUTOR_H
#define MPI3SNP_DISTRIBUTOR_H

#include <inttypes.h>
#include <vector>

class Distributor {
public:
    Distributor(unsigned int p_num, unsigned int p_id, unsigned long count);

    void Get_pairs(unsigned int t_num, unsigned int t_id,
                            std::vector<std::pair<uint32_t, uint32_t>> &values);

private:
    const unsigned long count;
    std::vector<uint32_t> distribution;
};

#endif //MPI3SNP_DISTRIBUTOR_H
