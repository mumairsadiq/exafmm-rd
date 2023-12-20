#include "kernel.h"
#include "mathfunc.h"

rtfmm::LaplaceKernel::LaplaceKernel()
{

}

void rtfmm::LaplaceKernel::p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, KernelType type)
{
    if(type == KernelType::naive)
    {
        printf("p2p_naive\n");
        for(int j = 0; j < cell_tar.body.number; j++)
        {
            Body3& btar = bs_tar[cell_tar.body.offset + j];
            real potential = 0;
            for(int i = 0; i < cell_src.body.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.body.offset + i];
                vec3r dx = btar.x - bsrc.x;
                real r = dx.r();
                real invr = r == 0 ? 0 : 1 / r;
                potential += bsrc.q * invr;
            }
            btar.p = potential;
        }
    }
    else if(type == KernelType::matrix)
    {
        printf("p2p matrix\n");
        int num_src = cell_src.body.number;
        int num_tar = cell_tar.body.number;

        std::vector<vec3r> x_src(num_src);
        for(int idx = 0; idx < num_src; idx++)
            x_src[idx] = bs_src[cell_src.body.offset + idx].x;

        std::vector<vec3r> x_tar(num_tar);
        for(int idx = 0; idx < num_tar; idx++)
            x_tar[idx] = bs_tar[cell_tar.body.offset + idx].x;

        std::vector<real> matrix_p2p(num_tar * num_src);
        get_p2p_matrix(x_src, x_tar, matrix_p2p);

        std::vector<real> q_src(num_src);
        for(int idx = 0; idx < num_src; idx++)
            q_src[idx] = bs_src[cell_src.body.offset + idx].q;

        std::vector<real> p_tar(num_tar);
        mat_vec(num_tar, num_src, matrix_p2p, q_src, p_tar, MathType::blas);

        for(int idx = 0; idx < num_tar; idx++)
            bs_tar[cell_tar.body.offset + idx].p = p_tar[idx];
    }
    else if(type == KernelType::fmm)
    {
        printf("p2p fmm\n");
    }
}

void rtfmm::LaplaceKernel::get_p2p_matrix(
    std::vector<vec3r>& x_src, 
    std::vector<vec3r>& x_tar, 
    std::vector<real>& matrix_p2p
)
{
    int num_src = x_src.size();
    int num_tar = x_tar.size();
    for(int j = 0; j < num_tar; j++)
    {
        vec3r xtar = x_tar[j];
        for(int i = 0; i < num_src; i++)
        {
            vec3r xsrc = x_src[i];
            vec3r dx = xtar - xsrc;
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            matrix_p2p[j * num_src + i] = invr;
        }
    }
}