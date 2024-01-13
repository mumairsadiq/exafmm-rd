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
Bodies3 generate_random_bodies(int num, real r, vec3r offset = vec3r(0,0,0), int seed = 0);

/**
 * @brief extract x from bodies
 * 
 * @param bs bodies
 * 
 * @return vector of x
 */
std::vector<vec3r> get_bodies_x(Bodies3& bs, Range range, vec3r offset = vec3r(0,0,0));

/**
 * @brief extract q from bodies
 * 
 * @param bs bodies
 * 
 * @return vector of q
 */
Matrix get_bodies_q(Bodies3& bs, Range range);


void set_boides_p(Bodies3& bs, Matrix& ps, Range range);

void set_boides_f(Bodies3& bs, Matriv& fs, Range range);

void add_boides_p(Bodies3& bs, Matrix& ps, Range range);

void add_boides_f(Bodies3& bs, Matriv& fs, Range range);


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
void print_bodies(const Bodies3& bs, int num = -1, int offset = 0, std::string name = "bodies");


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
BodyCompareResult compare(const Bodies3& bs1, const Bodies3& bs2, std::string name1, std::string name2, int num_compare = -1);

/**
 * @brief RT
 */
Bodies3 sort_bodies_by_idx(const Bodies3& bs);

}