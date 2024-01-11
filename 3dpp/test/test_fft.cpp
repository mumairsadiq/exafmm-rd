#include <fftw3.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "type.h"
#include "surface.h"
#include <map>

static std::map<int,rtfmm::vec3i> get_surface_conv_map(int p)
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


static std::vector<rtfmm::real> get_G_matrix(std::vector<rtfmm::vec3r>& grid, int N)
{
    std::vector<rtfmm::real> G;
    for(int k = 0; k < N; k++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int i = 0; i < N; i++)
            {
                rtfmm::real r = grid[k * N * N + j * N + i].r();
                rtfmm::real invr = r == 0 ? 0 : 1 / r;
                G.push_back(invr);
            }
        }
    }
    return G;
}

static std::vector<rtfmm::real> get_Q_matrix(std::vector<rtfmm::real> surface_q, int N, std::map<int,rtfmm::vec3i> surf_conv_map)
{
    std::vector<rtfmm::real> Q(N * N * N);
    int surface_number = surface_q.size();
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx = surf_conv_map[i];
        Q[idx[2] * N * N + idx[1] * N + idx[0]] = surface_q[i];
    }
    return Q;
}

static std::vector<rtfmm::real> naive(std::vector<rtfmm::real> qs, std::vector<rtfmm::vec3r> xs, rtfmm::vec3r offset)
{
    int N = qs.size();
    std::vector<rtfmm::real> ps(N);
    for(int j = 0; j < N; j++)
    {
        rtfmm::vec3r xj = xs[j];
        rtfmm::real p = 0;
        for(int i = 0; i < N; i++)
        {
            rtfmm::vec3r xi = xs[i] + offset;
            rtfmm::real r = (xi - xj).r();
            p += qs[i] / r;
        }
        ps[j] = p;
    }
    return ps;
}

static void print_real_grid(std::vector<rtfmm::real> grid, int N)
{
    for(int k = 0; k < N; k++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int i = 0; i < N; i++)
            {
                rtfmm::real x = grid[k * N * N + j * N + i];
                printf("%.4f  ",x);
            }
            printf("\n");
        }
        printf("\n\n");
    }
}

int main(int argc, char* argv[])
{
    int P = argc > 1 ? atoi(argv[1]) : 6;
    rtfmm::real r = 3.14 * 1.05;
    rtfmm::vec3r offset(1,2,3);
    int surface_number = rtfmm::get_surface_point_num(P);
    printf("P = %d, surf = %d, r = %.4f, offset = (%.2f,%.2f,%.2f)\n", P, surface_number, r, offset[0], offset[1], offset[2]);
    std::vector<rtfmm::vec3r> xs = rtfmm::get_surface_points(P, r, offset);
    std::vector<rtfmm::real> qs(surface_number);
    for(int i = 0; i < surface_number; i++)
    {
        qs[i] = drand48() - 0.5;
    }
    std::vector<rtfmm::real> ps_naive = naive(qs, xs, offset);

    rtfmm::real delta = 2 * r / (P - 1);
    int N = 2 * P - 1;
    rtfmm::real gsize = r * 2;
    std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, offset);
    std::vector<rtfmm::real> G = get_G_matrix(grid, N);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);
    std::vector<rtfmm::real> Q = get_Q_matrix(qs, N, surf_conv_map);

    std::vector<fftw_complex> Gk(N * N * N);
    std::vector<fftw_complex> Qk(N * N * N);
    std::vector<fftw_complex> pk(N * N * N);
    std::vector<rtfmm::real> p_fft_grid(N * N * N);
    fftw_plan plan_G = fftw_plan_dft_r2c_3d(N, N, N, G.data(), Gk.data(), FFTW_ESTIMATE);
    fftw_plan plan_Q = fftw_plan_dft_r2c_3d(N, N, N, Q.data(), Qk.data(), FFTW_ESTIMATE);
    fftw_execute(plan_G);
    fftw_execute(plan_Q);
    for(int i = 0; i < N * N * N; i++)
    {
        rtfmm::real a = Gk[i][0];
        rtfmm::real b = Gk[i][1];
        rtfmm::real c = Qk[i][0];
        rtfmm::real d = Qk[i][1];
        pk[i][0] = a * c - b * d;
        pk[i][1] = a * d + b * c;
    }
    fftw_plan plan_P = fftw_plan_dft_c2r_3d(N, N, N, pk.data(), p_fft_grid.data(), FFTW_ESTIMATE);
    fftw_execute(plan_P);
    for(int i = 0; i < N * N * N; i++)
    {
        p_fft_grid[i] /= N * N * N;
    }
    
    //std::vector<rtfmm::real> naive_grid = get_Q_matrix(ps_naive, N, surf_conv_map);
    //print_real_grid(naive_grid, N);
    //print_real_grid(p_fft_grid, N);

    std::vector<rtfmm::real> ps_fft(surface_number);
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx3 = surf_conv_map[i] + rtfmm::vec3i(P-1,P-1,P-1);
        ps_fft[i] = p_fft_grid[idx3[2] * N * N + idx3[1] * N + idx3[0]];
    }

    // compare
    rtfmm::real err = 0;
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::real p_naive = ps_naive[i];
        rtfmm::real p_fft = ps_fft[i];
        err += std::pow(p_naive - p_fft, 2);
        //printf("%.4f  %.4f\n", p_fft, p_naive);
    }
    err = std::sqrt(err / surface_number);
    std::cout<<"err = "<<err<<std::endl; 

    return 0;
}