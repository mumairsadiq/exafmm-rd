#include "kernel.h"
#include "mathfunc.h"
#include "surface.h"
#include "argument.h"
#include <omp.h>
#include <fftw3.h>

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

void rtfmm::LaplaceKernel::p2m_precompute(int P, Bodies3& bs_src, Cell3& cell_src)
{
    /* get source to check matrix */
    std::vector<vec3r> x_check = get_surface_points(P, cell_src.r * 2.95, cell_src.x);
    std::vector<vec3r> x_src = get_bodies_x(bs_src, cell_src.brange);
    Matrix s2c = get_p2p_matrix(x_src, x_check);

    /* get check potential */
    Matrix q_src = get_bodies_q(bs_src, cell_src.brange);
    Matrix p_check = mat_vec_mul(s2c, q_src);

    /* get equivalent charge */
    real scale = std::pow(0.5, cell_src.depth);
    Matrix UTb = mat_vec_mul(UT_p2m_precompute, p_check);
    Matrix SinvUTb = mat_vec_mul(Sinv_p2m_precompute, UTb);
    cell_src.q_equiv = mat_vec_mul(V_p2m_precompute, SinvUTb);
    mat_scale(cell_src.q_equiv, scale);
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

void rtfmm::LaplaceKernel::m2m_precompute(int P, Cell3& cell_parent, Cells3& cs)
{
    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];
        Matrix q_equiv_parent = mat_vec_mul(matrix_m2m[cell_child.octant], cell_child.q_equiv);
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


void rtfmm::LaplaceKernel::m2m_img_precompute(int P, Cell3& cell_parent, Cells3& cs, real cycle)
{
    Cell3& cell_child = cs[cell_parent.crange.offset];
    for(int octant = 0; octant < 27; octant++)
    {
        Matrix q_equiv_parent = mat_vec_mul(matrix_m2m_img[octant], cell_child.q_equiv);
        cell_parent.q_equiv = mat_mat_add(cell_parent.q_equiv, q_equiv_parent);
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


void rtfmm::LaplaceKernel::m2l_fft(int P, Cell3& cell_src, Cell3& cell_tar, vec3r offset)
{
    assert_exit(cell_src.depth == cell_tar.depth, "depth not equal");
    real r = cell_src.r * 1.05;
    vec3r relative_pos = cell_src.x + offset - cell_tar.x;
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, r, relative_pos);

    rtfmm::real delta = 2 * r / (P - 1);
    int N = 2 * P - 1;
    rtfmm::real gsize = r * 2;
    std::vector<rtfmm::vec3r> grid = get_conv_grid(N, gsize, delta, relative_pos);
    std::vector<rtfmm::real> G = get_G_matrix(grid, N);
    std::map<int,rtfmm::vec3i> surf_conv_map = get_surface_conv_map(P);
    std::vector<rtfmm::real> Q = get_Q_matrix(cell_src.q_equiv, N, surf_conv_map);

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

    int surface_number = rtfmm::get_surface_point_num(P);
    Matrix ps_fft(surface_number, 1);
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx3 = surf_conv_map[i] + rtfmm::vec3i(P-1,P-1,P-1);
        ps_fft.d[i] = p_fft_grid[idx3[2] * N * N + idx3[1] * N + idx3[0]];
    }

    fftw_destroy_plan(plan_P);
    fftw_destroy_plan(plan_G);
    fftw_destroy_plan(plan_Q);

    cell_tar.p_check = mat_mat_add(cell_tar.p_check, ps_fft);
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

void rtfmm::LaplaceKernel::l2l_precompute(int P, Cell3& cell_parent, Cells3& cs)
{
    for(int i = 0; i < cell_parent.crange.number; i++)
    {
        Cell3& cell_child = cs[cell_parent.crange.offset + i];
        Matrix p_check_child = mat_vec_mul(matrix_l2l[cell_child.octant], cell_parent.p_check);
        cell_child.p_check = mat_mat_add(cell_child.p_check, p_check_child);
    }
}

void rtfmm::LaplaceKernel::l2l_img_precompute(int P, Cell3& cell_parent, Cells3& cs)
{
    Cell3& cell_child = cs[cell_parent.crange.offset];
    Matrix p_check_child = mat_vec_mul(matrix_l2l_img, cell_parent.p_check);
    cell_child.p_check = mat_mat_add(cell_child.p_check, p_check_child);
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

void rtfmm::LaplaceKernel::l2p_precompute(int P, Bodies3& bs_tar, Cell3& cell_tar)
{
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_tar.r * 2.95, cell_tar.x);

    /* get equivalent charge */
    real scale = std::pow(0.5, cell_tar.depth);
    Matrix UTb = mat_vec_mul(UT_l2p_precompute, cell_tar.p_check);
    Matrix SinvUTb = mat_vec_mul(Sinv_l2p_precompute, UTb);
    Matrix q_equiv = mat_vec_mul(V_l2p_precompute, SinvUTb);
    mat_scale(q_equiv, scale);

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

void rtfmm::LaplaceKernel::precompute(int P, real r0)
{
    if(verbose) printf("precompute\n");

    /* p2m */
    std::vector<vec3r> x_check_up = get_surface_points(P, r0 * 2.95);
    std::vector<vec3r> x_equiv_up = get_surface_points(P, r0 * 1.05);
    Matrix e2c_up_precompute = get_p2p_matrix(x_equiv_up, x_check_up);

    Matrix U, S, VT;
    svd(e2c_up_precompute, U, S, VT);
    UT_p2m_precompute = transpose(U);
    V_p2m_precompute = transpose(VT);
    Sinv_p2m_precompute = pseudo_inverse(S);

    /* m2m */
    std::vector<vec3r>& x_equiv_up_parent = x_equiv_up;
    std::vector<vec3r>& x_check_up_parent = x_check_up;
    Matrix pe2pc_up_precompute = get_p2p_matrix(x_equiv_up_parent, x_check_up_parent);
    svd(pe2pc_up_precompute, U, S, VT);
    Matrix UT_m2m_precompute = transpose(U);
    Matrix V_m2m_precompute = transpose(VT);
    Matrix Sinv_m2m_precompute = pseudo_inverse(S);
    for(int octant = 0; octant < 8; octant++)
    {
        vec3r offset_child = Tree::get_child_cell_x(vec3r(0,0,0), r0, octant, 0);
        std::vector<vec3r> x_equiv_child_up = get_surface_points(P, r0 / 2 * 1.05, offset_child);
        Matrix ce2pc_up_precompute = get_p2p_matrix(x_equiv_child_up, x_check_up_parent);

        /* (VSinv)(UTG) */
        Matrix m1 = mat_mat_mul(UT_m2m_precompute, ce2pc_up_precompute);
        Matrix m2 = mat_mat_mul(V_m2m_precompute, Sinv_m2m_precompute);
        Matrix m3 = mat_mat_mul(m2, m1);
        matrix_m2m[octant] = m3;
    }

    /* m2m image */
    for(int octant = 0; octant < 27; octant++)
    {
        vec3r offset_child = Tree::get_child_cell_x(vec3r(0,0,0), r0, octant, 1);
        std::vector<vec3r> x_equiv_child_up = get_surface_points(P, r0 / 3 * 1.05, offset_child);
        Matrix ce2pc_up_precompute = get_p2p_matrix(x_equiv_child_up, x_check_up_parent);

        /* (VSinv)(UTG) */
        Matrix m1 = mat_mat_mul(UT_m2m_precompute, ce2pc_up_precompute);
        Matrix m2 = mat_mat_mul(V_m2m_precompute, Sinv_m2m_precompute);
        Matrix m3 = mat_mat_mul(m2, m1);
        matrix_m2m_img[octant] = m3;
    }

    /* l2l */
    std::vector<vec3r>& x_equiv_down = x_check_up;
    std::vector<vec3r>& x_check_down = x_equiv_up;
    Matrix pe2pc_down_precompute = get_p2p_matrix(x_equiv_down, x_check_down);
    svd(pe2pc_down_precompute, U, S, VT);
    Matrix UT_l2l_precompute = transpose(U);
    Matrix V_l2l_precompute = transpose(VT);
    Matrix Sinv_l2l_precompute = pseudo_inverse(S);

    std::vector<vec3r>& x_equiv_parent_down = x_equiv_down;
    for(int octant = 0; octant < 8; octant++)
    {
        vec3r offset_child = Tree::get_child_cell_x(vec3r(0,0,0), r0, octant, 0);
        std::vector<vec3r> x_equiv_child_down = get_surface_points(P, r0 / 2 * 1.05, offset_child);
        Matrix pe2cc_down_precompute = get_p2p_matrix(x_equiv_parent_down, x_equiv_child_down);

        /* (G(VSinv))UT */
        Matrix m1 = mat_mat_mul(V_l2l_precompute, Sinv_l2l_precompute);
        Matrix m2 = mat_mat_mul(pe2cc_down_precompute, m1);
        Matrix m3 = mat_mat_mul(m2, UT_l2l_precompute);

        /* G((VSinv)UT)*/
        /*Matrix m1 = mat_mat_mul(V_l2l_precompute, Sinv_l2l_precompute);
        Matrix m2 = mat_mat_mul(m1, UT_l2l_precompute);
        Matrix m3 = mat_mat_mul(pe2cc_down_precompute, m2);*/

        matrix_l2l[octant] = m3;
    }

    /* l2l image */
    {
        vec3r offset_child = vec3r(0,0,0);
        std::vector<vec3r> x_equiv_child_down = get_surface_points(P, r0 / 3 * 1.05, offset_child);
        Matrix pe2cc_down_precompute = get_p2p_matrix(x_equiv_parent_down, x_equiv_child_down);

        /* (G(VSinv))UT */
        Matrix m1 = mat_mat_mul(V_l2l_precompute, Sinv_l2l_precompute);
        Matrix m2 = mat_mat_mul(pe2cc_down_precompute, m1);
        Matrix m3 = mat_mat_mul(m2, UT_l2l_precompute);

        matrix_l2l_img = m3;
    }

    /* l2p */
    Matrix e2c_down_precompute = get_p2p_matrix(x_equiv_down, x_check_down);

    svd(e2c_down_precompute, U, S, VT);
    UT_l2p_precompute = transpose(U);
    V_l2p_precompute = transpose(VT);
    Sinv_l2p_precompute = pseudo_inverse(S);
}

std::vector<rtfmm::real> rtfmm::LaplaceKernel::get_G_matrix(std::vector<rtfmm::vec3r>& grid, int N)
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

std::vector<rtfmm::real> rtfmm::LaplaceKernel::get_Q_matrix(Matrix& surface_q, int N, std::map<int,rtfmm::vec3i> surf_conv_map)
{
    std::vector<rtfmm::real> Q(N * N * N);
    int surface_number = surface_q.m * surface_q.n;
    for(int i = 0; i < surface_number; i++)
    {
        rtfmm::vec3i idx = surf_conv_map[i];
        Q[idx[2] * N * N + idx[1] * N + idx[0]] = surface_q[i];
    }
    return Q;
}

std::map<int,rtfmm::vec3i> rtfmm::LaplaceKernel::get_surface_conv_map(int p)
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