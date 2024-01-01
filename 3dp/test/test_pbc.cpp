#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"
#include "ewald.h"
#include <omp.h>

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3dp_test_pbc"<<std::endl;

    rtfmm::Argument args(argc, argv);
    args.show();

    omp_set_num_threads(args.th_num);
    printf("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x, args.seed);

    /* solve by FMM */
    rtfmm::LaplaceFMM fmm(bs, args);
    rtfmm::Bodies3 res_fmm = fmm.solve();

    /* solve by ewald */
    rtfmm::EwaldSolver ewald(bs, args);
    rtfmm::Bodies3 res_ewald = ewald.solve();

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    rtfmm::LaplaceKernel kernel;
    rtfmm::Cell3 cell_src;
    cell_src.brange = {0, args.n};
    rtfmm::Cell3 cell_tar;
    cell_tar.brange = {0, args.num_compare};
    TIME_BEGIN(direct);
    int dm = (std::pow(3, args.images) - 1) / 2;
    printf("dm = %d\n", dm);
    for(int pz = -dm; pz <= dm; pz++)
    {
        printf("pz = %d\n", pz);
        for(int py = -dm; py <= dm; py++)
        {
            for(int px = -dm; px <= dm; px++)
            {
                kernel.p2p(bs, res_direct, cell_src, cell_tar, rtfmm::vec3r(px,py,pz) * args.cycle);
            }   
        }
    }
    rtfmm::dipole_correction(res_direct, args.cycle);
    if(args.timing) {TIME_END(direct);}

    /* compare */
    if(rtfmm::verbose)
    {
        rtfmm::print_bodies(res_fmm, 3, 0, "fmm");
        rtfmm::print_bodies(res_direct, 3, 0, "direct");
        rtfmm::print_bodies(res_ewald, 3, 0, "ewald");
    }
    rtfmm::compare(res_fmm, res_direct, "FMM", "Direct", args.num_compare).show();
    rtfmm::compare(res_fmm, res_ewald, "FMM", "Ewald", args.num_compare).show();
    rtfmm::compare(res_direct, res_ewald, "Direct", "Ewald", args.num_compare).show();

    return 0;
}