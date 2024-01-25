#include "surface.h"

std::vector<rtfmm::vec3r> rtfmm::get_surface_points(int p, real r, vec3r x, int dir)
{
    int num = get_surface_point_num(p);
    std::vector<rtfmm::vec3r> points;
    if(dir == 0)
    {
        for(int i = 0; i < p; i++)
        {
            for(int j = 0; j < p; j++)
            {
                for(int k = 0; k < p; k++)
                {
                    if(i == 0 || i == p - 1 || j == 0 || j == p - 1 || k == 0 || k == p - 1)
                    {
                        points.push_back(rtfmm::vec3r(-1.0 + i * 2.0 / (p - 1),-1.0 + j * 2.0 / (p - 1),-1.0 + k * 2.0 / (p - 1)));
                    }
                }
            }
        }
    }
    else
    {
        for(int k = 0; k < p; k++)
        {
            for(int j = 0; j < p; j++)
            {
                for(int i = 0; i < p; i++)
                {
                    if(i == 0 || i == p - 1 || j == 0 || j == p - 1 || k == 0 || k == p - 1)
                    {
                        points.push_back(rtfmm::vec3r(-1.0 + i * 2.0 / (p - 1),-1.0 + j * 2.0 / (p - 1),-1.0 + k * 2.0 / (p - 1)));
                    }
                }
            }
        }
    }
    for(int i = 0; i < points.size(); i++)
    {
        vec3r& p = points[i];
        p = p * r + x;
    }
    assert_exit(points.size() == num, "surface point number error");

    return points;
}

std::vector<rtfmm::vec3r> rtfmm::get_conv_grid(int grid_len, rtfmm::real gmin, rtfmm::real delta, rtfmm::vec3r offset)
{
    std::vector<rtfmm::vec3r> grid;
    for(int k = 0; k < grid_len; k++)
    {
        for(int j = 0; j < grid_len; j++)
        {
            for(int i = 0; i < grid_len; i++)
            {
                rtfmm::real x = -gmin + i * delta;
                rtfmm::real y = -gmin + j * delta;
                rtfmm::real z = -gmin + k * delta;
                rtfmm::vec3r point(x,y,z);
                grid.push_back(point - offset);
            }
        }
    }
    return grid;
}