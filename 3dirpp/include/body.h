#pragma once
#include "type.h"
#include <vector>
#include <algorithm>

namespace rtfmm
{

struct Body
{
    int idx;
    real q;
    vec3r x;

    real p;
    vec3r f;

    friend std::ostream& operator<<(std::ostream& os, const Body& b)
    {
        os << "<Body> "
        << "idx=" << b.idx
        << ", q=" << b.q
        << ", x=" << b.x
        << ", p=" << b.p
        << ", f=" << b.f;
        return os;
    }
};

using Bodies = std::vector<Body>;

/**
 * @brief Generate random bodies in 3d box.
 * 
 * @param num number of bodies
 * @param r half size of cubic box
 * @param offset box center
 * 
 * @return bodies with zero net charge 
 */
Bodies generate_random_bodies(int num, real r, vec3r offset = vec3r(0,0,0), int seed = 0, int zero_netcharge = 1);

inline Bodies sort_bodies_by_idx(const Bodies& bs)
{
    Bodies bodies = bs;
    std::sort(
        bodies.begin(), 
        bodies.end(),
        [](const Body& a, const Body& b)
        {
            return a.idx < b.idx; 
        }
    );

    return bodies;
}

}