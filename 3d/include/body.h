#pragma once
#include "type.h"
#include <vector>

namespace rtfmm
{

struct Body3
{
    int idx;
    real q;
    real p;
    vec3r x;
    vec3r f;
};

using Bodies3 = std::vector<Body3>;

/**
 * @brief Generate random bodies in 3d box.
 * 
 * @param num number of bodies
 * @param r half size of cubic box
 * @param offset box center
 * 
 * @return bodies with zero net charge 
 */
Bodies3 generate_random_bodies(int num, real r, vec3r offset = vec3r(0,0,0));

/**
 * @brief extract x from bodies
 * 
 * @param bs bodies
 * 
 * @return vector of x
 */
std::vector<vec3r> get_bodies_x(Bodies3& bs, Range range);

/**
 * @brief extract q from bodies
 * 
 * @param bs bodies
 * 
 * @return vector of q
 */
std::vector<real> get_bodies_q(Bodies3& bs, Range range);


/**
 * @brief set potential of bodies from ps
 * 
 * @param bs bodies
 * @param ps potentials
 */
void set_boides_p(Bodies3& bs, std::vector<real>& ps, Range range);

/**
 * @brief set force of bodies from ps
 * 
 * @param bs bodies
 * @param fs potentials
 */
void set_boides_f(Bodies3& bs, std::vector<vec3r>& fs, Range range);


/**
 * @brief Print one body.
 * 
 * @param b body
 */
void print_body(const Body3& b);


/**
 * @brief Print many bodies.
 * 
 * @param bs bodies
 * @param num if -1, print all; else, print min(num,bs.size()) bodies
 * @param offset offset
 */
void print_bodies(const Bodies3& bs, int num = -1, int offset = 0);

}