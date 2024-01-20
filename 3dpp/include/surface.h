#pragma once
#include "type.h"
#include <vector>

namespace rtfmm
{

int get_surface_point_num(int p);

/**
*@brief generate surface points warpping a box
*@param p number of points in each direction
*@param r half of surface cube size
*@param x center of box
*/
std::vector<vec3r> get_surface_points(int p, real r = 1, vec3r x = vec3r(0,0,0), int dir = 0);

std::vector<rtfmm::vec3r> get_conv_grid(int grid_len, rtfmm::real gmin, rtfmm::real delta, rtfmm::vec3r offset);

}