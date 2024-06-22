#include "type.h"
#include "argument.h"
#include "surface.h"
#include "mathfunc.h"
#include "kernel.h"
#include <functional>

int main(int argc, char* argv[])
{
    using namespace rtfmm;
    title("rtfmm_3dpp_test_matmul");

    Argument args(argc, argv);
    args.show();

    std::vector<vec3r> x_check_up = get_surface_points(args.P, args.r * 2.95, 1);
    std::vector<vec3r> x_equiv_up = get_surface_points(args.P, args.r * 1.05, 0);
    LaplaceKernel kernel;
    Matrix e2c = kernel.get_p2p_matrix(x_equiv_up, x_check_up);
    /*for(int j = 0; j < e2c.m; j++)
    {
        for(int i = 0; i < e2c.n; i++)
        {
            e2c.d[j * e2c.n + i] = drand48();
        }
    }*/
    //print_matrix(e2c);
    Matrix U, S, VT;
    svd(e2c, U, S, VT);
    
    Matrix UT = transpose(U);
    Matrix V = transpose(VT);
    Matrix Sinv = pseudo_inverse(S);
    //print_matrix(S);
    Bodies3 bs = generate_random_bodies(args.n, args.r, args.x, 5, args.zero_netcharge);
    Cell3 cell_src;
    cell_src.idx = 0;
    cell_src.depth = 0;
    cell_src.x = args.x;
    cell_src.r = args.r;
    cell_src.brange = Range(0, args.n);
    std::vector<vec3r> x_src = get_bodies_x(bs, cell_src.brange);
    Matrix s2c = kernel.get_p2p_matrix(x_src, x_check_up);
    Matrix q_src = get_bodies_q(bs, cell_src.brange);
    Matrix p_check = mat_vec_mul(s2c, q_src);
    real scale = std::pow(0.5, cell_src.depth);

    if(0)
    {
        Matrix ones = identity(e2c.m, e2c.n);
        Matrix VSinv = mat_mat_mul(V, Sinv);
        Matrix VSinvUT = mat_mat_mul(VSinv, UT);
        Matrix i1 = mat_mat_mul(VSinvUT, e2c);
        std::cout<<"error1 = "<<matrix_L2(i1, ones)<<std::endl;

        Matrix SinvUT = mat_mat_mul(Sinv, UT);
        VSinvUT = mat_mat_mul(V, SinvUT);
        Matrix i2 = mat_mat_mul(VSinvUT, e2c);
        std::cout<<"error2 = "<<matrix_L2(i2, ones)<<std::endl;
    }

    if(0)
    {
        /* (VSinv)(UTb) */
        Matrix UTb = mat_vec_mul(UT, p_check);
        Matrix VSinv = mat_mat_mul(V, Sinv);
        Matrix res1 = mat_vec_mul(VSinv, UTb, scale);
        std::cout<<res1.d[0]<<std::endl;

        /* ((VSinv)UT)b) */
        Matrix VSinvUT = mat_mat_mul(VSinv, UT);
        Matrix res2 = mat_vec_mul(VSinvUT, p_check, scale);
        std::cout<<res2.d[0]<<std::endl;

        /* V(Sinv(UTb) */
        Matrix SinvUTb = mat_mat_mul(Sinv, UTb);
        Matrix res3 = mat_mat_mul(V, SinvUTb);
        std::cout<<res3.d[0]<<std::endl;

        /* V((SinvUT)b) */
        Matrix SinvUT = mat_mat_mul(Sinv, UT);
        SinvUTb = mat_mat_mul(SinvUT, p_check);
        Matrix res4 = mat_mat_mul(V, SinvUTb);
        std::cout<<res4.d[0]<<std::endl;

        /* (V(SinvUT))b) */
        SinvUT = mat_mat_mul(Sinv, UT);
        VSinvUT = mat_mat_mul(V, SinvUT);
        Matrix res5 = mat_vec_mul(VSinvUT, p_check, scale);
        std::cout<<res5.d[0]<<std::endl;

        std::cout<<"L2 = "<<matrix_L2(res1, res2)<<std::endl;
        std::cout<<"L2 = "<<matrix_L2(res1, res3)<<std::endl;
        std::cout<<"L2 = "<<matrix_L2(res1, res4)<<std::endl;
        std::cout<<"L2 = "<<matrix_L2(res1, res5)<<std::endl;
    }

    if(1)
    {
        auto mmm = mat_mat_mul;
        //std::function<Matrix(const Matrix&,const Matrix&)> mmm = mat_mat_mul_naive;
        /* A((VSinv)(UTb)) = b */
        Matrix UTb = mmm(UT, p_check);
        Matrix VSinv = mmm(V, Sinv);
        Matrix VSinvUTb = mmm(VSinv, UTb);
        Matrix res1 = mmm(e2c, VSinvUTb);
        std::cout<<"error1 = "<<matrix_L2(res1, p_check)<<std::endl;

        /* A((VSinv)UT)b) = b */
        Matrix VSinvUT = mmm(VSinv, UT);
        VSinvUTb = mmm(VSinvUT, p_check);
        Matrix res2 = mmm(e2c, VSinvUTb);
        std::cout<<"error2 = "<<matrix_L2(res2, p_check)<<std::endl;
    }

    return 0;
}