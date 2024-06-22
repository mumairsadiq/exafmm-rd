#pragma once
#include "type.h"
#include <vector>
#include <map>
#include <ostream>

namespace rtfmm
{

inline int get_surface_point_num(int p)
{
    return 6 * (p - 1) * (p - 1) + 2;
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

/**
 * @brief
 * N = 2 * P - 1;
 * N3 = N * N * N;
 * N_freq = N * N * (N / 2 + 1);
 * delta = 2 * r / (P - 1) * 1.05;
 * gsize = r * 2 * 1.05;
*/
struct conv_grid_setting
{
    conv_grid_setting(int P, real r)
    {
        N = 2 * P - 1; // check P.41 and P.42 of kernel_independent.pptx
        N3 = N * N * N;
        N_freq = N * N * (N / 2 + 1);
        delta = 2 * r / (P - 1) * 1.05;
        gsize = r * 2 * 1.05;
    }
    int N; // number of points per dim
    int N3; // total number of points(as the name implies, N * N * N)
    /**
     * @brief total number of points in Fourier space, N * N * (N / 2 + 1)
     * @note according to fft's property, 
     * when inputs are real numbers, the out[N/2,N] is conj to out[1,N/2], 
     * therefore FFTW left out[N/2,N] zero, 
     * so when use r2c or c2r, we also store half of full frequences
    */
    int N_freq;
    real delta; // interval distance of grid points per dim
    real gsize; // size of grid per dim

    friend std::ostream &operator<<(std::ostream& os, conv_grid_setting& conv_grid)
    {
        os << "conv_grid(N = " << conv_grid.N << ", N3 = "<<conv_grid.N3<<", N_freq = "<<conv_grid.N_freq<<", delta = "<<conv_grid.delta<<", gsize = "<<conv_grid.gsize<<")";
        return os;
    }
};

}