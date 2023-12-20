#pragma once
#include "type.h"
#include "body.h"
#include "cell.h"

namespace rtfmm
{

class LaplaceFMM
{
public:
    LaplaceFMM();

    LaplaceFMM(
        Bodies3& bs,
        int P_,
        vec3r x_,
        real r_
    );

    void solve();

private:
    Bodies3 bodies;
    int num_body;
    int P;
    vec3r x;
    real r;
    Cells3 cells;
};

}