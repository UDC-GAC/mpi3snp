//
// Created by christian on 23/01/17.
//

#ifndef BOOLVECTOR_H
#define BOOLVECTOR_H

#include <cstdlib>
#include <cstring>

class BoolVector{
public:
    BoolVector();
    ~BoolVector();
    size_t trueCount();
    size_t falseCount();
    void push_back(bool value);
    bool get(size_t pos);
    size_t size();
    void resize(size_t newSize);
    inline bool &operator [](size_t pos){
        return array[pos];
    }
private:
    static const int INITIAL_SIZE = 1024;
    bool *array;
    size_t arraySize, count, tCount, fCount;
};

#endif //BOOLVECTOR_H
