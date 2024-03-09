#include "type.h"
#include "tree.h"
#include "argument.h"
#include "body.h"
#include "tree.h"
#include "traverser.h"
#include <omp.h>
#include "kernel.h"

int main(int argc, char* argv[])
{
    rtfmm::title("rtfmm/3dirpp/test/test_p2p");
    rtfmm::Argument args(argc, argv);
    args.show();

    omp_set_dynamic(0);
    omp_set_num_threads(args.th_num);
    printf("# of threads = %d\n", omp_get_max_threads());

    /* prepare bodies */
    rtfmm::Bodies bs1 = rtfmm::generate_random_bodies(args.n, 1, rtfmm::vec3r(0,0,0), args.seed, args.zero_netcharge);
    rtfmm::Bodies bs2 = rtfmm::generate_random_bodies(args.n, 1, rtfmm::vec3r(0,0,0), args.seed, args.zero_netcharge);

    rtfmm::Cell ctar;
    ctar.brange = rtfmm::vec2i(0, args.n);
    for(int i = 0; i < bs1.size(); i++)
    {
        ctar.bs.push_back(bs1[i]);
    }
    rtfmm::Cell csrc;
    csrc.brange = rtfmm::vec2i(0, args.n);
    for(int i = 0; i < bs2.size(); i++)
    {
        csrc.bs.push_back(bs2[i]);
    }
    rtfmm::Cells cs;
    cs.push_back(ctar);
    cs.push_back(csrc);
    std::vector<int> p2ps;
    p2ps.push_back(1);
    rtfmm::LaplaceKernel kernel;
    kernel.p2p(cs, p2ps, ctar);
    rtfmm::dipole_correction(ctar.bs, args.cycle);

    kernel.direct(bs2,bs1,args.images,args.cycle);
    rtfmm::dipole_correction(bs1,args.cycle);

    if(rtfmm::verbose)
    {
        if(args.enable_fmm) rtfmm::print_bodies(ctar.bs, args.print_body_number, 0, "fmm");
        if(args.enable_direct) rtfmm::print_bodies(bs1, args.print_body_number, 0, "direct");
    }
    rtfmm::compare(bs1, ctar.bs, "direct", "p2p").show();
    
    std::cout<<"test pass!"<<std::endl;

    return 0;
}