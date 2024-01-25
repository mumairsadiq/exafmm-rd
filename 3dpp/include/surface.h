#pragma once
#include "type.h"
#include <vector>
#include <map>

namespace rtfmm
{

inline int get_surface_point_num(int p)
{
    return 6 * (p - 1) * (p - 1) + 2;
}

/**
 * @brief 2 * P - 1
 * @note check P.41 and P.42 of kernel_independent.pptx
*/
inline int get_conv_grid_point_num_per_dim(int P)
{
    int N = 2 * P - 1;
    return N;
}

/**
 * @brief N * N * (N / 2 + 1), where N = 2 * P - 1
 * @note according to fft's properity, 
 * when inputs are real numbers, the out[N/2,N] is conj to out[1,N/2], 
 * therefore FFTW left out[N/2,N] zero, 
 * so when use r2c or c2r, we also store half of full frequences
*/
inline int get_nfreq(int P)
{
    int N = 2 * P - 1;
    int N_freq = N * N * (N / 2 + 1); 
    return N_freq;
}

/**
*@brief generate surface points warpping a box
*@param p number of points in each direction
*@param r half of surface cube size
*@param x center of box
*/
std::vector<vec3r> get_surface_points(int p, real r = 1, vec3r x = vec3r(0,0,0), int dir = 0);

std::vector<rtfmm::vec3r> get_conv_grid(int grid_len, rtfmm::real gmin, rtfmm::real delta, rtfmm::vec3r offset);

/**
 * @brief index-of-surface-point -> coordinate-in-conv-grid
*/
inline std::map<int,rtfmm::vec3i> get_surface_conv_map(int p)
{
    int num = rtfmm::get_surface_point_num(p);
    std::map<int,rtfmm::vec3i> map;
    int cnt = 0;
    for(int i = 0; i < p; i++)
    {
        for(int j = 0; j < p; j++)
        {
            for(int k = 0; k < p; k++)
            {
                if(i == 0 || i == p - 1 || j == 0 || j == p - 1 || k == 0 || k == p - 1)
                {
                    map[cnt] = rtfmm::vec3i(i,j,k);
                    cnt++;
                }
            }
        }
    }
    return map;
}

}