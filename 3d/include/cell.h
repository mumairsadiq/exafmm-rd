#pragma once
#include "type.h"
#include "body.h"

namespace rtfmm
{

struct Cell3
{
    int idx;
    int depth;
    real r;
    vec3r x;
    OffsetAndNumber child;
    OffsetAndNumber body;

    matrix q_equiv;
    matrix p_check;
};

using Cells3 = std::vector<Cell3>;

}