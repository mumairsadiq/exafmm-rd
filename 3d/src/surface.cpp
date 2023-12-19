#include "surface.h"

std::vector<rtfmm::vec3r> rtfmm::get_surface_points(int p)
{
    std::vector<rtfmm::vec3r> points;
    for(int i = 0; i < p; i++)
    {
        for(int j = 0; j < p; j++)
        {
            for(int k = 0; k < p; k++)
            {
                if(i == 0 || i == p - 1 || j == 0 || j == p - 1 || k == 0 || k == p - 1)
                {
                    points.push_back(rtfmm::vec3r(-1 + i * 2 / (p - 1),-1 + j * 2 / (p - 1),-1 + k * 2 / (p - 1)));
                }
            }
        }
    }
    return points;
}