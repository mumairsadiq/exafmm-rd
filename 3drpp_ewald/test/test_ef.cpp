#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include "ewald.h"
#include "direct.h"
#include <omp.h>
#include "timer.h"
#include "epkifmm.h"

int main(int argc, char* argv[])
{
    rtfmm::title(argv[0]);

    rtfmm::Argument args(argc, argv);
    args.show(args.res_filepath);

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    if(rtfmm::verbose) printf("# of threads = %d\n", omp_get_max_threads());

    rtfmm::Timer timer;
    timer.begin("TOTAL");

    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::Bodies3 res_fmm;
    if(args.enable_fmm)
    {
        rtfmm::EpkiFMM epkifmm(bs, args);
        TIME_BEGIN(FMM);
        timer.begin("EPKIFMM");
        res_fmm = epkifmm.solve();
        if(args.timing) {TIME_END(FMM);timer.end("EPKIFMM");}
    }

    /* solve by ewald */
    rtfmm::Bodies3 res_ewald;
    if(args.enable_ewald)
    {
        rtfmm::EwaldSolver ewald(bs, bs, args);
        TIME_BEGIN(EWALD);
        timer.begin("EWALD");
        res_ewald = ewald.solve();
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_ewald);
        if(args.timing) {TIME_END(EWALD);timer.end("EWALD");}
    }

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    if(args.enable_direct)
    {
        rtfmm::LaplaceKernel kernel;
        rtfmm::Cell3 cell_src;
        cell_src.brange = {0, args.n};
        rtfmm::Cell3 cell_tar;
        cell_tar.brange = {0, args.num_compare};
        TIME_BEGIN(DIRECT);
        timer.begin("DIRECT");
        kernel.direct(res_direct, res_direct, args.images, args.cycle);
        if(args.dipole_correction)
            rtfmm::dipole_correction(res_direct, args.cycle);
        if(args.divide_4pi)
            rtfmm::scale_bodies(res_direct);
        if(args.timing) {TIME_END(DIRECT);timer.end("DIRECT");}
    }

    timer.end("TOTAL");
    timer.save(args.res_filepath);

    /* compare */
    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "epkifmm");
        if(args.enable_direct) rtfmm::print_bodies(res_direct, args.print_body_number, 0, "direct");
        if(args.enable_ewald) rtfmm::print_bodies(res_ewald, args.print_body_number, 0, "ewald");
    }
    if(args.enable_fmm && args.enable_direct)
        rtfmm::compare(res_fmm, res_direct, "EPKIFMM", "Direct", args.num_compare).show(args.res_filepath);
    if(args.enable_fmm && args.enable_ewald)
        rtfmm::compare(res_fmm, res_ewald, "EPKIFMM", "Ewald", args.num_compare).show(args.res_filepath);
    if(args.enable_direct && args.enable_ewald)
        rtfmm::compare(res_direct, res_ewald, "Direct", "Ewald", args.num_compare).show(args.res_filepath);

    return 0;
}