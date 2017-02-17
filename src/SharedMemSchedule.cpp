//
// Created by christian on 14/02/17.
//

#include "SharedMemSchedule.h"

SharedMemSchedule::SharedMemSchedule(size_t size, int numProc) :
        size(size), numProc(numProc) {
    int i, j;
    int blockSize = size / numProc,
            remainder = size % numProc;

    // Allocate space
    blockDiv.reserve(numProc);
    comb.reserve(numProc * (numProc + 1) / 2);
    proc.reserve(numProc * (numProc + 1) / 2);

    // Create basic division in SNP contiguous blocks
    for (i = 0; i < numProc; i++) {
        int beginning = i > 0 ? blockDiv[i - 1]->x + blockDiv[i - 1]->xlen : 0;
        int size = blockSize + (remainder-- > 0 ? 1 : 0);
        blockDiv.push_back(new Block(beginning, size));
    }

    // Create block distribution
    for (i = 0; i < numProc; i++) {
        for (j = i; j < numProc; j++) {
            std::vector<Block *> item;
            item.push_back(blockDiv[i]);
            item.push_back(blockDiv[j]);
            comb.push_back(item);
        }
    }

    // Assign different blocks to each process
    for (i = 0, j = 0; i < comb.size(); i++) {
        proc.push_back(j++);
        j %= numProc;
    }
}

SharedMemSchedule::~SharedMemSchedule() {
    // Deallocate Block objects
    for (auto it = blockDiv.begin(); it != blockDiv.end(); it++){
        delete *it;
    }
}

std::vector<std::vector<Block *>> SharedMemSchedule::getBlocks(int procId) {
    std::vector<std::vector<Block *>> res;
    for (int i = 0; i < comb.size(); i++) {
        if (proc[i] == procId) {
            res.push_back(comb[i]);
        }
    }

    return res;
}

int SharedMemSchedule::whoHas(std::vector<Block *> comb) {
    // Everyone has every block
    return -1;
}