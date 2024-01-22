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
    int seed;
    int num_compare;
    int ewald_ksize;
    int th_num;

    int dipole_correction;
    int zero_netcharge;
    int print_body_number;
    int divide_4pi;

    int enable_fmm;
    int enable_direct;
    int enable_ewald;
    int check_tree;
    int use_fft;
    int use_precompute;
    int use_simd;
};

}