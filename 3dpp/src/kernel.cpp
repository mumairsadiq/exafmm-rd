#include "kernel.h"
#include "mathfunc.h"
#include "surface.h"
#include "argument.h"
#include <omp.h>

rtfmm::LaplaceKernel::LaplaceKernel()
{
    
}

void rtfmm::LaplaceKernel::p2p(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    #pragma omp parallel for
    for(int j = 0; j < cell_tar.brange.number; j++)
    {
        Body3& btar = bs_tar[cell_tar.brange.offset + j];
        real potential = 0;
        vec3r force(0,0,0);
        for(int i = 0; i < cell_src.brange.number; i++)
        {
            Body3& bsrc = bs_src[cell_src.brange.offset + i];
            vec3r dx = btar.x - bsrc.x - offset;
            real r = dx.r();
            real invr = r == 0 ? 0 : 1 / r;
            potential += bsrc.q * invr;
            force += bsrc.q * invr * invr * invr * (-dx);
        }
        btar.p += potential;
        btar.f += force;
    }
}

void rtfmm::LaplaceKernel::p2p_matrix(Bodies3& bs_src, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange, offset);
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    Matrix matrix_p2p = get_p2p_matrix(x_src, x_tar);
    Matrix p_tar = mat_vec_mul(matrix_p2p, q_src);
    Matriv f_tar = get_force_naive(x_src, x_tar, q_src);
    add_boides_p(bs_tar, p_tar, cell_tar.brange);
    add_boides_f(bs_tar, f_tar, cell_tar.brange);
}

void rtfmm::LaplaceKernel::p2m(int P, Bodies3& bs_src, Cell3& cell_src)
{
    /* get source to check matrix */
    std::vector<vec3r> x_check = get_surface_points(P, cell_src.r * 2.95, cell_src.x);
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange);
    Matrix s2c = get_p2p_matrix(x_src, x_check);

    /* get check potential */
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    Matrix p_check = mat_vec_mul(s2c, q_src);

    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x);
    Matrix e2c = get_p2p_matrix(x_equiv, x_check);

    /* get equivalent charge */
    cell_src.q_equiv = linear_equation_system_svd(e2c, p_check);
}


void rtfmm::LaplaceKernel::m2m(int P, Cell3& cell_parent, Cells3& cs)
{
    /* parent equiv to parent check matrix */
    std::vector<rtfmm::vec3r> x_equiv_parent = get_surface_points(P, cell_parent.r * 1.05, cell_parent.x);
    std::vector<rtfmm::vec3r> x_check_parent = get_surface_points(P, cell_parent.r * 2.95, cell_parent.x);
    Matrix pe2pc = get_p2p_matrix(x_equiv_parent, x_check_parent);

    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];

        /* child equivalent to parent check */
        std::vector<rtfmm::vec3r> x_equiv_child = get_surface_points(P, cell_child.r * 1.05, cell_child.x);
        Matrix ce2pc = get_p2p_matrix(x_equiv_child, x_check_parent);
        Matrix p_check_parent = mat_vec_mul(ce2pc, cell_child.q_equiv);

        /* parent check potential to parent equivalent charge */
        Matrix q_equiv_parent = linear_equation_system_svd(pe2pc, p_check_parent);
        cell_parent.q_equiv = mat_mat_add(cell_parent.q_equiv, q_equiv_parent);
    }
}

void rtfmm::LaplaceKernel::m2m_img(int P, Cell3& cell_parent, Cells3& cs, real cycle)
{
    /* parent equiv to parent check matrix */
    std::vector<rtfmm::vec3r> x_equiv_parent = get_surface_points(P, cell_parent.r * 1.05, cell_parent.x);
    std::vector<rtfmm::vec3r> x_check_parent = get_surface_points(P, cell_parent.r * 2.95, cell_parent.x);
    Matrix pe2pc = get_p2p_matrix(x_equiv_parent, x_check_parent);

    for(int pz = -1; pz <= 1; pz++)
    {
        for(int py = -1; py <= 1; py++)
        {
            for(int px = -1; px <= 1; px++)
            {
                Cell3& cell_child = cs[cell_parent.crange.offset];

                /* child equivalent to parent check */
                std::vector<rtfmm::vec3r> x_equiv_child = get_surface_points(P, cell_child.r * 1.05, cell_child.x + vec3r(px,py,pz) * cycle);
                Matrix ce2pc = get_p2p_matrix(x_equiv_child, x_check_parent);
                Matrix p_check_parent = mat_vec_mul(ce2pc, cell_child.q_equiv);

                /* parent check potential to parent equivalent charge */
                Matrix q_equiv_parent = linear_equation_system_svd(pe2pc, p_check_parent);
                cell_parent.q_equiv = mat_mat_add(cell_parent.q_equiv, q_equiv_parent);
            }
        }
    }
}


void rtfmm::LaplaceKernel::m2l(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    /* get src equivalent to tar check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x + offset);
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix e2c = get_p2p_matrix(x_equiv, x_check);

    // source cell equivalent charge to target cell check potential
    Matrix p_check = mat_vec_mul(e2c, cell_src.q_equiv);
    cell_tar.p_check = mat_mat_add(cell_tar.p_check, p_check);
}


void rtfmm::LaplaceKernel::l2l(int P, Cell3& cell_parent, Cells3& cs)
{
    /* parent check to parent equiv */
    std::vector<rtfmm::vec3r> x_equiv_parent = get_surface_points(P, cell_parent.r * 2.95, cell_parent.x);
    std::vector<rtfmm::vec3r> x_check_parent = get_surface_points(P, cell_parent.r * 1.05, cell_parent.x);
    Matrix pe2pc = get_p2p_matrix(x_equiv_parent, x_check_parent);
    Matrix q_equiv_parent = linear_equation_system_svd(pe2pc, cell_parent.p_check);

    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];

        /* parent equiv to child check */
        std::vector<rtfmm::vec3r> x_check_child = get_surface_points(P, cell_child.r * 1.05, cell_child.x);
        Matrix pe2cc = get_p2p_matrix(x_equiv_parent, x_check_child);
        Matrix p_check_child = mat_vec_mul(pe2cc, q_equiv_parent);
        cell_child.p_check = mat_mat_add(cell_child.p_check, p_check_child);
    }
}


void rtfmm::LaplaceKernel::l2p(int P, Bodies3& bs_tar, Cell3& cell_tar)
{
    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_tar.r * 2.95, cell_tar.x);
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix e2c = get_p2p_matrix(x_equiv, x_check);

    /* get equivalent charge */
    Matrix q_equiv = linear_equation_system_svd(e2c, cell_tar.p_check);

    /* get equiv to tar matrix */
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix e2t = get_p2p_matrix(x_equiv, x_tar);

    /* get target potential and force */
    Matrix p_tar = mat_vec_mul(e2t, q_equiv);
    Matriv f_tar = get_force_naive(x_equiv, x_tar, q_equiv);
    add_boides_p(bs_tar, p_tar, cell_tar.brange);
    add_boides_f(bs_tar, f_tar, cell_tar.brange);
}

void rtfmm::LaplaceKernel::m2p(int P, Bodies3& bs_tar, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    std::vector<rtfmm::vec3r> x_equiv_src = get_surface_points(P, cell_src.r * 1.05, cell_src.x + offset);
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar, cell_tar.brange);
    Matrix e2t = get_p2p_matrix(x_equiv_src, x_tar);
    Matrix p_tar = mat_vec_mul(e2t, cell_src.q_equiv);
    Matriv f_tar = get_force_naive(x_equiv_src, x_tar, cell_src.q_equiv);
    add_boides_p(bs_tar, p_tar, cell_tar.brange);
    add_boides_f(bs_tar, f_tar, cell_tar.brange);
}


void rtfmm::LaplaceKernel::p2l(int P, Bodies3& bs_src, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange, offset);
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    std::vector<rtfmm::vec3r> x_check_tar = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    Matrix s2c = get_p2p_matrix(x_src, x_check_tar);
    Matrix p_check = mat_vec_mul(s2c, q_src);
    cell_tar.p_check = mat_mat_add(cell_tar.p_check, p_check);
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

void rtfmm::dipole_correction(Bodies3& bs, real cycle)
{
    if(verbose) std::cout<<"dipole correction"<<std::endl;
    int num_body = bs.size();
    real coef = 4 * M_PI / (3 * cycle * cycle * cycle);
    vec3r dipole(0,0,0);
    for (int i = 0; i < num_body; i++) 
    {
        dipole += bs[i].x * bs[i].q;
    }
    real dnorm = dipole.norm();
    for (int i = 0; i < num_body; i++) 
    { 
        bs[i].p -= coef * dnorm / num_body / bs[i].q; 
        bs[i].f -= coef * dipole;
    }
}

void rtfmm::LaplaceKernel::precompute()
{
    
}