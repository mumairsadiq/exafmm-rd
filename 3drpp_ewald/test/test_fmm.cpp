#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include <omp.h>
#include "timer.h"

int main(int argc, char* argv[])
{
    rtfmm::title(argv[0]);

    rtfmm::Argument args(argc, argv);
    args.show(args.res_filepath);

    omp_set_num_threads(args.th_num);
    if(rtfmm::verbose) printf("# of threads = %d\n", omp_get_max_threads());

    rtfmm::Timer timer;

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed, args.zero_netcharge);

    /* solve by FMM */
    rtfmm::LaplaceFMM fmm(bs, args);
    TIME_BEGIN(FMM);
    timer.begin("FMM");
    rtfmm::Bodies3 res_fmm = fmm.solve();
    if(args.timing) {TIME_END(FMM);timer.end("FMM");}

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    rtfmm::LaplaceKernel kernel;
    rtfmm::Cell3 cell;
    cell.brange = {0, args.n};
    TIME_BEGIN(direct);
    timer.begin("direct");
    kernel.p2p(bs, res_direct, cell, cell, rtfmm::vec3r(0,0,0), args.use_simd);
    if(args.divide_4pi)
        rtfmm::scale_bodies(res_direct);
    if(args.dipole_correction)
        rtfmm::dipole_correction(res_direct, args.cycle);
    if(args.timing) {TIME_END(direct);timer.end("direct");}

    timer.save(args.res_filepath);

    /* compare */
    if(rtfmm::verbose)
    {
        rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "res_fmm");
        rtfmm::print_bodies(res_direct, args.print_body_number, 0, "res_direct");
    }
    rtfmm::compare(res_fmm, res_direct, "FMM", "Direct").show(args.res_filepath);

    return 0;
}