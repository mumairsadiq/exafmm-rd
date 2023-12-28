#include "type.h"
#include "body.h"
#include "tree.h"
#include "fmm.h"
#include "argument.h"

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3d_test_fmm"<<std::endl;

    rtfmm::Argument args(argc, argv);
    args.show();

    /* prepare bodies */
    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(args.n, args.r, args.x);

    /* solve by FMM */
    rtfmm::LaplaceFMM fmm(bs, args);
    rtfmm::Bodies3 res_fmm = fmm.solve();

    /* solve directly */
    rtfmm::Bodies3 res_direct = bs;
    rtfmm::LaplaceKernel kernel;
    rtfmm::Cell3 cell;
    cell.brange = {0, args.n};
    kernel.p2p(bs, res_direct, cell, cell);

    /* compare */
    //rtfmm::print_bodies(res_fmm, -1, 0, "res_fmm");
    //rtfmm::print_bodies(res_direct, -1, 0, "res_direct");
    rtfmm::compare(res_fmm, res_direct, "FMM", "Direct").show();

    return 0;
}