#include "kernel.h"
#include "mathfunc.h"
#include "surface.h"

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
            vec3r force(0,0,0);
            for(int i = 0; i < cell_src.body.number; i++)
            {
                Body3& bsrc = bs_src[cell_src.body.offset + i];
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
    else if(type == KernelType::matrix)
    {
        printf("p2p matrix\n");
        int num_src = cell_src.body.number;
        int num_tar = cell_tar.body.number;

        std::vector<vec3r> x_src = get_bodies_x(bs_src);
        std::vector<vec3r> x_tar = get_bodies_x(bs_tar);
        std::vector<real> q_src = get_bodies_q(bs_src);
        std::vector<real> matrix_p2p(num_tar * num_src);
        std::vector<real> p_tar(num_tar);

        get_p2p_matrix(x_src, x_tar, matrix_p2p);
        mat_vec(num_tar, num_src, matrix_p2p, q_src, p_tar, MathType::blas);
        set_boides_p(bs_tar, p_tar);
    }
    else if(type == KernelType::fmm)
    {
        printf("p2p fmm\n");
    }
}


void rtfmm::LaplaceKernel::p2m(int P, Bodies3& bs_src, Cell3& cell_src)
{
    /* get source to check matrix */
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_src.r * 2.95, cell_src.x);
    int num_check = x_check.size();
    int num_src = cell_src.body.number;
    std::vector<vec3r> x_src = get_bodies_x(bs_src);
    matrix matrix_s2c(num_check * num_src);
    get_p2p_matrix(x_src, x_check, matrix_s2c);

    /* get check potential */
    matrix q_src = get_bodies_q(bs_src);
    matrix p_check(num_check * 1);
    mat_vec(num_check, num_src, matrix_s2c, q_src, p_check, MathType::blas);

    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x);
    int num_equiv = x_equiv.size();
    matrix matrix_e2c(num_check * num_equiv);
    get_p2p_matrix(x_equiv, x_check, matrix_e2c);

    /* get equivalent charge */
    int num_s = std::min(num_check, num_equiv);
    matrix matrix_e2c_U(num_check * num_check);
    matrix matrix_e2c_S(num_s * 1);
    matrix matrix_e2c_VT(num_equiv * num_equiv);
    svd(num_equiv, num_check, matrix_e2c, matrix_e2c_U, matrix_e2c_S, matrix_e2c_VT);
    matrix UT = transpose(num_check, num_check, matrix_e2c_U);
    matrix V = transpose(num_equiv, num_equiv, matrix_e2c_VT);
    matrix Sinv(num_check * num_equiv);
    for(int j = 0; j < num_check; j++)
    {
        for(int i = 0; i < num_equiv; i++)
        {
            if(i == j && j < num_s) Sinv[j * num_equiv + i] = 1.0 / matrix_e2c_S[j];
            else Sinv[j * num_equiv + i] = 0;
        }
    }
    Sinv = transpose(num_check, num_equiv, Sinv);
    matrix VSinvUTp(num_equiv * 1); // q_equiv = V * S^-1 * U^T * p_check
    matrix UTp(num_check * 1);
    matrix SinvUTp(num_equiv * 1);
    mat_vec(num_check, num_check, UT, p_check, UTp, MathType::blas);
    mat_vec(num_equiv, num_check, Sinv, UTp, SinvUTp, MathType::blas);
    mat_vec(num_equiv, num_check, V, SinvUTp, VSinvUTp, MathType::blas);

    cell_src.q_equiv = VSinvUTp;

    /* check result */
    /*matrix p_check2(num_check * 1);
    mat_vec(num_check, num_equiv, matrix_e2c, cell_src.q_equiv, p_check2, MathType::blas);
    real err = 0;
    for(int i = 0; i < num_check; i++)
    {
        err += std::pow(p_check[i] - p_check2[i], 2);
    }
    err = std::sqrt(err / num_check);
    std::cout<<"error = "<<err<<std::endl;*/
}


void rtfmm::LaplaceKernel::m2l(int P, Cell3& cell_src, Cell3& cell_tar)
{
    /* get src equivalent to tar check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_src.r * 1.05, cell_src.x);
    int num_equiv = x_equiv.size();
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    int num_check = x_check.size();
    matrix matrix_e2c(num_check * num_equiv);
    get_p2p_matrix(x_equiv, x_check, matrix_e2c);

    // source cell equivalent charge to target cell potential
    mat_vec(num_check, num_equiv, matrix_e2c, cell_src.q_equiv, cell_tar.p_check, MathType::blas);
}


void rtfmm::LaplaceKernel::l2p(int P, Bodies3& bs_tar, Cell3& cell_tar)
{
    /* get equivalent to check matrix */
    std::vector<rtfmm::vec3r> x_equiv = get_surface_points(P, cell_tar.r * 2.95, cell_tar.x);
    int num_equiv = x_equiv.size();
    std::vector<rtfmm::vec3r> x_check = get_surface_points(P, cell_tar.r * 1.05, cell_tar.x);
    int num_check = x_check.size();
    matrix matrix_e2c(num_check * num_equiv);
    get_p2p_matrix(x_equiv, x_check, matrix_e2c);

    /* get equivalent charge */
    int num_s = std::min(num_check, num_equiv);
    matrix matrix_e2c_U(num_check * num_check);
    matrix matrix_e2c_S(num_s * 1);
    matrix matrix_e2c_VT(num_equiv * num_equiv);
    svd(num_equiv, num_check, matrix_e2c, matrix_e2c_U, matrix_e2c_S, matrix_e2c_VT);
    matrix UT = transpose(num_check, num_check, matrix_e2c_U);
    matrix V = transpose(num_equiv, num_equiv, matrix_e2c_VT);
    matrix Sinv(num_check * num_equiv);
    for(int j = 0; j < num_check; j++)
    {
        for(int i = 0; i < num_equiv; i++)
        {
            if(i == j && j < num_s) Sinv[j * num_equiv + i] = 1.0 / matrix_e2c_S[j];
            else Sinv[j * num_equiv + i] = 0;
        }
    }
    Sinv = transpose(num_check, num_equiv, Sinv);
    matrix q_equiv(num_equiv * 1); // q_equiv = V * S^-1 * U^T * p_check
    matrix UTp(num_check * 1);
    matrix SinvUTp(num_equiv * 1);
    mat_vec(num_check, num_check, UT, cell_tar.p_check, UTp, MathType::blas);
    mat_vec(num_equiv, num_check, Sinv, UTp, SinvUTp, MathType::blas);
    mat_vec(num_equiv, num_check, V, SinvUTp, q_equiv, MathType::blas);

    /* get equiv to tar matrix */
    std::vector<vec3r> x_tar = get_bodies_x(bs_tar);
    int num_tar = x_tar.size();
    matrix matrix_e2t(num_tar * num_equiv);
    get_p2p_matrix(x_equiv, x_tar, matrix_e2t);

    /* get target potential */
    matrix p_tar(num_tar * 1);
    mat_vec(num_tar, num_equiv, matrix_e2t, q_equiv, p_tar, MathType::blas);
    set_boides_p(bs_tar, p_tar);
    std::vector<vec3r> f_tar(num_tar * 1);
    get_force_naive(x_equiv, x_tar, q_equiv, f_tar);
    set_boides_f(bs_tar, f_tar);
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

void rtfmm::LaplaceKernel::get_force_naive(
    std::vector<vec3r> x_src, 
    std::vector<vec3r> x_tar, 
    std::vector<real> q_src,
    std::vector<vec3r>& f_tar
)
{
    std::vector<vec3r> res(x_tar.size());
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
    f_tar = res;
}
