#include "type.h"
#include "body.h"
#include "cell.h"
#include "kernel.h"

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3d_test_m2l"<<std::endl;

    int p = argc > 1 ? atoi(argv[1]) : 4;
    int s = argc > 2 ? atoi(argv[2]) : 10;
    printf("[p = %d, s = %d]\n", p, s);

    rtfmm::real r = 1.0;
    int num_body_src = s;
    int num_body_tar = 5;
    rtfmm::vec3r x_src(-3,0,0);
    rtfmm::vec3r x_tar(3,0,0);

    rtfmm::Bodies3 bs_src = rtfmm::generate_random_bodies(num_body_src, r, x_src);
    rtfmm::Cell3 cell_src;
    cell_src.idx = 0;
    cell_src.depth = 0;
    cell_src.r = r;
    cell_src.x = x_src;
    cell_src.child = {0,0};
    cell_src.body = {0,num_body_src};

    rtfmm::Bodies3 bs_tar = rtfmm::generate_random_bodies(num_body_tar, r, x_tar);
    rtfmm::Cell3 cell_tar;
    cell_tar.idx = 0;
    cell_tar.depth = 0;
    cell_tar.r = r;
    cell_tar.x = x_tar;
    cell_tar.child = {0,0};
    cell_tar.body = {0,num_body_tar};
    
    rtfmm::LaplaceKernel kernel;

    /* FMM */
    kernel.p2m(p, bs_src, cell_src);
    kernel.m2l(p, cell_src, cell_tar);
    kernel.l2p(p, bs_tar, cell_tar);

    /* naive */
    rtfmm::Bodies3 bs_tar2 = bs_tar;
    rtfmm::Cell3 cell_tar2 = cell_tar;
    kernel.p2p(bs_src, bs_tar2, cell_src, cell_tar2, rtfmm::LaplaceKernel::KernelType::naive);

    /* compare */
    rtfmm::print_bodies(bs_tar, num_body_tar);
    rtfmm::print_bodies(bs_tar2, num_body_tar);
    rtfmm::real perr = 0.0, ferr = 0.0f;
    rtfmm::real fdif = 0.0, fnrm = 0.0;
    for(int i = 0; i < num_body_tar; i++)
    {
        perr += std::pow(bs_tar[i].p - bs_tar2[i].p, 2);
        fdif += (bs_tar[i].f - bs_tar2[i].f).norm();
        fnrm += bs_tar2[i].f.norm();
    }
    perr = std::sqrt(perr / num_body_tar);
    ferr = std::sqrt(fdif / fnrm);
    std::cout<<"L2p = "<<perr<<", L2f = "<<ferr<<std::endl;

    return 0;
}