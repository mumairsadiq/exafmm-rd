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
std::vector<vec3r> get_surface_points(int p, real r = 1, vec3r x = vec3r(0,0,0));

}