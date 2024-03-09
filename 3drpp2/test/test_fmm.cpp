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

    /* when i is zero, simply enlarge the fmm's radius */
    if(args.images == 0) args.r += args.rega;

    /* solve by FMM */
    rtfmm::LaplaceFMM fmm(bs, args);
    TIME_BEGIN(FMM);
    rtfmm::Bodies3 res_fmm = fmm.solve();
    if(args.timing) {TIME_END(FMM);}

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    rtfmm::LaplaceKernel kernel;
    rtfmm::Cell3 cell_src;
    cell_src.brange = {0, args.n};
    rtfmm::Cell3 cell_tar;
    cell_tar.brange = {0, args.num_compare};
    TIME_BEGIN(DIRECT);
    kernel.direct(res_direct, res_direct, args.images, args.cycle);
    if(args.dipole_correction)
        rtfmm::dipole_correction(res_direct, args.cycle);
    if(args.divide_4pi)
        rtfmm::scale_bodies(res_direct);
    if(args.timing) {TIME_END(DIRECT);}

    /* compare */
    if(rtfmm::verbose)
    {
        rtfmm::print_bodies(res_fmm, args.print_body_number, 0, "res_fmm");
        rtfmm::print_bodies(res_direct, args.print_body_number, 0, "res_direct");
    }
    rtfmm::compare(res_fmm, res_direct, "FMM", "Direct").show();

    return 0;
}