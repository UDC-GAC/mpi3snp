//
// Created by christian on 06/02/17.
//

#ifndef BLOCK_H
#define BLOCK_H

#include <stdlib.h>

class Block {
public:
    const size_t x, xlen;

    inline Block(size_t x, size_t xlen) :
            x(x), xlen(xlen) {
    }

    friend bool operator==(const Block &lhs, const Block &rhs) {
        return lhs.x == rhs.x && lhs.xlen == rhs.xlen;
    }

    friend bool operator!=(const Block &lhs, const Block &rhs) {
        return !(lhs.x == rhs.x && lhs.xlen == rhs.xlen);
    }

    friend bool operator<=(const Block &lhs, const Block &rhs) {
        return lhs.x <= rhs.x;
    }

    friend bool operator<(const Block &lhs, const Block &rhs) {
        return lhs.x < rhs.x;
    }

    friend bool operator>(const Block &lhs, const Block &rhs) {
        return lhs.x > rhs.x;
    }

    friend bool operator>=(const Block &lhs, const Block &rhs) {
        return lhs.x >= rhs.x;
    }
};



#endif //BLOCK_H
