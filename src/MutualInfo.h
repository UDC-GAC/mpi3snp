/*
 * MutualInfo.h
 *
 *  Created on: 19/09/2014
 *      Author: jorge
 */

#ifndef MUTUALINFO_H_
#define MUTUALINFO_H_

#include <inttypes.h>
#include <string>

/*
 * Structure for the outputs with the result of the mutual informations and the SNPs associated
 */

struct MutualInfo {
    MutualInfo() {
        _id1 = 0;
        _id2 = 0;
        _id3 = 0;
        _mutualInfoValue = 0.0;
    }

    uint32_t _id1;
    uint32_t _id2;
    uint32_t _id3;
    float _mutualInfoValue;

    bool operator<(const MutualInfo mi) const { return _mutualInfoValue < mi._mutualInfoValue; }

    inline std::string To_string() {
        return std::to_string(_id1) + " " + std::to_string(_id2) + " " + std::to_string(_id3) + " " +
               std::to_string(_mutualInfoValue);
    }
};

#endif // MUTUALINFO_H_
