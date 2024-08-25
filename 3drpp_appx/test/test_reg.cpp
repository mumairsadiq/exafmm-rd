#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include "ewald.h"
#include "direct.h"
#include <omp.h>


int main(int argc, char* argv[])
{
    rtfmm::title("rtfmm_3dpp_test_reg");

    rtfmm::Argument args(argc, argv);
    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    if(rtfmm::verbose) printf("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r-0.05, args.x, args.seed, args.zero_netcharge);
    //rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r - args.rega * 2, args.x, args.seed, args.zero_netcharge);

    if(args.body0_idx != -1)
    {
        RTLOG("move body %d\n", args.body0_idx);
        bs[args.body0_idx].x = rtfmm::vec3r(args.x0, args.y0, args.z0);
        //bs[args.body0_idx].q = 0;
    }

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::LaplaceFMM fmm(bs, args);
        TIME_BEGIN(FMM);
        res_fmm = fmm.solve();
        if(args.timing) {TIME_END(FMM);}
    }

    /* solve by ewald */
    rtfmm::Bodies3 res_ewald;
    if(args.enable_ewald)
    {
        rtfmm::EwaldSolver ewald(bs, args);
        TIME_BEGIN(EWALD);
        res_ewald = ewald.solve();
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_ewald);
        if(args.timing) {TIME_END(EWALD);}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct;
    if(args.enable_direct)
    {
        rtfmm::Cell3 cell_tar;
        cell_tar.bodies = bs;
        rtfmm::Cell3 cell_src;
        rtfmm::make_source_cell(cell_src, bs);

        TIME_BEGIN(DIRECT);
        rtfmm::LaplaceKernel kernel;
        kernel.direct(cell_src, cell_tar, args.images, args.cycle);
        res_direct = cell_tar.bodies;
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END(DIRECT);}
    }

    /* compare */
    {
        if(args.enable_fmm)    rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
        if(args.enable_ewald)  rtfmm::print_bodies(res_ewald, args.print_body_number, 0, "ewald");
    }

    if(args.body0_idx != -1)
    {
        RTLOG("check : ");
        rtfmm::Body3& cb = res_fmm[args.check_body_idx];
        /*std::cout << "idx=" << cb.idx << ","
                   << "p=" << cb.p << ","
                   << "f=" << cb.f << std::endl;*/
        printf("idx = %d, p = %.12f\n", cb.idx, cb.p);
    }

    if(args.enable_fmm && args.enable_direct)
        rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare).show();
    if(args.enable_fmm && args.enable_ewald)
        rtfmm::compare(res_fmm, res_ewald, "FMM", "Ewald", args.num_compare).show();
    if(args.enable_direct && args.enable_ewald)
        rtfmm::compare(res_direct, res_ewald, "Direct", "Ewald", args.num_compare).show();

    return 0;
}