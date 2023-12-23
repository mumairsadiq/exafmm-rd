#pragma once

#include "type.h"

namespace rtfmm
{

class Argument
{
public:
    Argument(int argc, char* argv[])
    {
        r = 1.0;
        x = vec3r(0,0,0);
        P = argc > 1 ? atoi(argv[1]) : 4;
        n = argc > 2 ? atoi(argv[2]) : 1000;
    }
    void show()
    {
        printf("[P = %d, n = %d, r = %.4f, x = (%.3f,%.3f,%.3f)]\n", P, n, r, x[0], x[1], x[2]);
    }

public:
    vec3r x;
    real r;
    int n;
    int P;
};

}