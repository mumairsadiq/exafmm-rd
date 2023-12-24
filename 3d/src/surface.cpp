#include "surface.h"

int rtfmm::get_surface_point_num(int p)
{
    return 6 * (p - 1) * (p - 1) + 2;
}

std::vector<rtfmm::vec3r> rtfmm::get_surface_points(int p, real r, vec3r x)
{
    int num = get_surface_point_num(p);
    std::vector<rtfmm::vec3r> points;
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
    for(int i = 0; i < points.size(); i++)
    {
        vec3r& p = points[i];
        p = p * r + x;
    }
    if(points.size() != num)
    {
        printf("point number error ! %ld, %d\n", points.size(), num);
        exit(0);
    }
    return points;
}