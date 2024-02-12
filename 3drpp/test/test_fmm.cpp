#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include <omp.h>

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3dp_test_fmm"<<std::endl;

    rtfmm::Argument args(argc, argv);
    args.show();

    omp_set_num_threads(args.th_num);
    if(rtfmm::verbose) printf("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, 5, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::LaplaceFMM fmm(bs, args);
    TIME_BEGIN(FMM);
    rtfmm::Bodies3 res_fmm = fmm.solve();
    if(args.timing) {TIME_END(FMM);}

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    rtfmm::LaplaceKernel kernel;
    rtfmm::Cell3 cell;
    cell.brange = {0, args.n};
    cell.bs.clear();
    for(int i = 0; i < args.n; i++)
    {
        cell.bs.push_back(bs[i]);
        cell.ws.push_back(1.0);
    }
    TIME_BEGIN(direct);
    kernel.p2p(cell, cell, rtfmm::vec3r(0,0,0), args.use_simd);
    for(int i = 0; i < args.n; i++)
    {
        res_direct[i] = cell.bs[i];
    }
    /*std::vector<rtfmm::vec3r> xs = rtfmm::get_bodies_x(res_direct, cell.brange);
    std::vector<rtfmm::real> qs = rtfmm::get_bodies_q(res_direct, cell.brange).d;
    std::vector<rtfmm::real> ps(args.n);
    std::vector<rtfmm::vec3r> fs(args.n);
    rtfmm::LaplaceKernel::p2p_pf(xs,qs,xs,ps,fs);
    rtfmm::Matrix psmat(args.n, 1);
    psmat.d = ps;
    rtfmm::Matriv fsmat(args.n, 1);
    fsmat.d = fs;
    rtfmm::add_boides_p(res_direct,psmat,cell.brange);
    rtfmm::add_boides_f(res_direct,fsmat,cell.brange);*/
    if(args.divide_4pi)
        rtfmm::scale_bodies(res_direct);
    if(args.dipole_correction)
        rtfmm::dipole_correction(res_direct, args.cycle);
    if(args.timing) {TIME_END(direct);}

    /* compare */
    if(rtfmm::verbose)
    {
        rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "res_fmm");
        rtfmm::print_bodies(res_direct, args.print_body_number, 0, "res_direct");
    }
    rtfmm::compare(res_fmm, res_direct, "FMM", "Direct").show();

    return 0;
}