#pragma once

#include "type.h"
#include "cmdline.h"

namespace rtfmm
{

extern int verbose;

class Argument
{
public:
    Argument(int argc, char* argv[]);
    void show();

public:
    vec3r x;
    real r;
    int n;
    int P;
    int ncrit;
    int timing;
    int images;
    real cycle;
};

}