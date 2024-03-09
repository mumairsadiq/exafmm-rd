#pragma once
#include "type.h"
#include <vector>
#include <algorithm>
#include "math.h"

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

inline void print_bodies(const Bodies& bs, int num, int offset, std::string name)
{
    std::cout<<name<<":"<<std::endl;
    if(num == -1) num = bs.size();
    int s = std::min(num, (int)bs.size());
    for(int i = offset; i < offset + s; i++)
    {
        std::cout<<bs[i]<<std::endl;
    }
    printf("\n");
}

void dipole_correction(Bodies& bs, real cycle);

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

struct BodyCompareResult
{
    int num_compared;
    real rmsp;
    real rmsf;
    real l2p;
    real l2f;
    real epot1;
    real epot2;
    real l2e;
    std::string name1;
    std::string name2;
    void show();
};

/**
 * @brief compare two Bodies's potential and force to get information(rms-error, l2-error, potential-energy)
 * 
 * @param bs1 bodies1
 * @param bs2 bodies2
 * @param name1 name of bs1
 * @param name2 name of bs2
 * @param num_compare number of comparison(-1 means all)
 * 
 * @return result of body data comparison
 * 
 * @warning bs2 is used to calculate norm
 */
BodyCompareResult compare(const Bodies& bs1, const Bodies& bs2, std::string name1, std::string name2, int num_compare = -1);

std::vector<vec3r> get_bodies_x(Bodies& bs, vec2i range, vec3r offset);

Matrix get_bodies_q(Bodies& bs, vec2i range);

Matriv get_force_naive(
    std::vector<vec3r>& x_src, 
    std::vector<vec3r>& x_tar, 
    Matrix& q_src
);

}