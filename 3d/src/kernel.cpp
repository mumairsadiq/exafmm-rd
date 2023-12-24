#include "kernel.h"
#include "mathfunc.h"
#include "surface.h"

rtfmm::LaplaceKernel::LaplaceKernel()
{

}

void rtfmm::LaplaceKernel::p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar)
{
    for(int j = 0; j < cell_tar.brange.number; j++)
    {
        Body3& btar = bs_tar[cell_tar.brange.offset + j];
        real potential = 0;
        vec3r force(0,0,0);
        for(int i = 0; i < cell_src.brange.number; i++)
        {
            Body3& bsrc = bs_src[cell_src.brange.offset + i];
            vec3r dx = btar.x - bsrc.x;
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            potential += bsrc.q * invr;
            force += bsrc.q * invr * invr * invr * (-dx);
        }
        btar.p = potential;
        btar.f = force;
    }
}

void rtfmm::LaplaceKernel::p2p_matrix(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar)
{
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange);
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    Matrix matrix_p2p = get_p2p_matrix(x_src, x_tar);
    Matrix p_tar = mat_vec(matrix_p2p, q_src);
    Matriv f_tar = get_force_naive(x_src, x_tar, q_src);
    set_boides_p(bs_tar, p_tar, cell_tar.brange);
    set_boides_f(bs_tar, f_tar, cell_tar.brange);
}

void rtfmm::LaplaceKernel::p2m(int P, Bodies3& bs_src, Cell3& cell_src)
{
    /* get source to check matrix */
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_src.r * 2.95, cell_src.x);
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange);
    Matrix matrix_s2c = get_p2p_matrix(x_src, x_check);

    /* get check potential */
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    Matrix p_check = mat_vec(matrix_s2c, q_src);

    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x);
    Matrix matrix_e2c = get_p2p_matrix(x_equiv, x_check);

    /* get equivalent charge */
    cell_src.q_equiv = linear_equation_system_svd(matrix_e2c, p_check);
}


void rtfmm::LaplaceKernel::m2m(int P, Cell3& cell_src)
{

}


void rtfmm::LaplaceKernel::m2l(int P, Cell3& cell_src, Cell3& cell_tar)
{
    /* get src equivalent to tar check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x);
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix matrix_e2c = get_p2p_matrix(x_equiv, x_check);

    // source cell equivalent charge to target cell check potential
    cell_tar.p_check = mat_vec(matrix_e2c, cell_src.q_equiv);
}


void rtfmm::LaplaceKernel::l2l(int P, Cell3& cell_tar)
{

}


void rtfmm::LaplaceKernel::l2p(int P, Bodies3& bs_tar, Cell3& cell_tar)
{
    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_tar.r * 2.95, cell_tar.x);
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix matrix_e2c = get_p2p_matrix(x_equiv, x_check);

    /* get equivalent charge */
    Matrix q_equiv = linear_equation_system_svd(matrix_e2c, cell_tar.p_check);

    /* get equiv to tar matrix */
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix matrix_e2t = get_p2p_matrix(x_equiv, x_tar);

    /* get target potential and force */
    Matrix p_tar = mat_vec(matrix_e2t, q_equiv);
    Matriv f_tar = get_force_naive(x_equiv, x_tar, q_equiv);
    set_boides_p(bs_tar, p_tar, cell_tar.brange);
    set_boides_f(bs_tar, f_tar, cell_tar.brange);
}


rtfmm::Matrix rtfmm::LaplaceKernel::get_p2p_matrix(
    std::vector<vec3r>& x_src, 
    std::vector<vec3r>& x_tar
)
{
    int num_src = x_src.size();
    int num_tar = x_tar.size();
    Matrix matrix_p2p(num_tar, num_src);
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
    return matrix_p2p;
}

rtfmm::Matriv rtfmm::LaplaceKernel::get_force_naive(
    std::vector<vec3r> x_src, 
    std::vector<vec3r> x_tar, 
    Matrix q_src
)
{
    Matriv res(x_tar.size(), 1);
    for(int j = 0; j < x_tar.size(); j++)
    {   
        vec3r force(0,0,0);
        for(int i = 0; i < x_src.size(); i++)
        {
            vec3r dx = x_tar[j] - x_src[i];
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            force += q_src[i] * invr * invr * invr * (-dx);
        }
        res[j] = force;
    }
    return res;
}
