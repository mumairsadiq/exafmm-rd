#pragma once
#include "type.h"
#include <math.h>
#include <vector>

namespace rtfmm
{

template<typename T>
std::vector<T> exclusive_scan(std::vector<T>& d)
{
    std::vector<T> res(d.size());
    T sum = 0;
    size_t ss = d.size();
    for(size_t i = 0; i < ss; i++)
    {
        res[i] = sum;
        sum += d[i];
    }
    return res;
}

}