//
// Created by christian on 06/02/17.
//

#ifndef BLOCK_H
#define BLOCK_H

#include <stdlib.h>

class Block {
public:
    const size_t x, xlen;

    inline Block(size_t x, size_t xlen):
            x(x), xlen(xlen)
    {
    }
};


#endif //BLOCK_H
