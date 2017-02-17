//
// Created by christian on 23/01/17.
//

#include "BoolVector.h"

BoolVector::BoolVector(){
    array = NULL;
    arraySize = 0;
    count = 0;
    tCount = 0;
    fCount = 0;
}

BoolVector::~BoolVector(){
    delete array;
}

size_t BoolVector::falseCount() {
    return fCount;
}

size_t BoolVector::trueCount() {
    return tCount;
}

void BoolVector::push_back(bool value){
    if (count == arraySize){
        reserve(arraySize * 2 + 1);
    }

    array[count++] = value;
    if (value){
        tCount++;
    } else {
        fCount++;
    }
}

bool BoolVector::get(size_t pos){
    return array[pos];
}

size_t BoolVector::size() {
    return count;
}

void BoolVector::reserve(size_t newSize){
    bool *temp = new bool[newSize];
    memcpy(temp, array, sizeof(bool) * arraySize);
    delete array;
    array = temp;
    arraySize = newSize;
}