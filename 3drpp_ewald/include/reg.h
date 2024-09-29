#pragma once
#include "type.h"

namespace rtfmm
{
namespace reg
{
inline real reg_w(real x)
{
    return 0.25 * (2 + 3 * x - x * x * x);
}

inline real get_w_single(real dx, real R, real rega)
{
    real r = std::abs(dx) - R + rega;
    real x;
    if(r <= 0) x = 1;
    else if(r > 2 * rega) x = -1;
    else x = 1 - r / rega;
    return reg_w(x);
}

inline real get_w(vec3r dx, real R, real rega)
{
    real w = 1;
    for(int d = 0; d < 3; d++)
    {
        w *= get_w_single(dx[d], R, rega);
    }
    return w;
}

inline vec3r get_w_xyz(vec3r dx, real R, real rega)
{
    rtfmm::vec3r res;
    for(int d = 0; d < 3; d++)
    {
        res[d] = get_w_single(dx[d], R, rega);
    }
    return res;
}
}
}