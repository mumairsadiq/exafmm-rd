#include "type.h"
#include "body.h"
#include "cell.h"
#include "kernel.h"

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3d_test_p2p"<<std::endl;

    rtfmm::real r = 1.0;
    int num_body_src = 1000;
    int num_body_tar = 5;
    rtfmm::vec3r x_src(-3,0,0);
    rtfmm::vec3r x_tar(3,0,0);

    rtfmm::Bodies3 bs_src = rtfmm::generate_random_bodies(num_body_src, r, x_src);
    rtfmm::Cell3 cell_src;
    cell_src.idx = 0;
    cell_src.depth = 0;
    cell_src.r = r;
    cell_src.x = {0,0,0};
    cell_src.child = {0,0};
    cell_src.body = {0,num_body_src};

    rtfmm::Bodies3 bs_tar = rtfmm::generate_random_bodies(num_body_tar, r, x_tar);
    rtfmm::Cell3 cell_tar;
    cell_tar.idx = 0;
    cell_tar.depth = 0;
    cell_tar.r = r;
    cell_tar.x = {0,0,0};
    cell_tar.child = {0,0};
    cell_tar.body = {0,num_body_tar};
    
    rtfmm::LaplaceKernel kernel;
    kernel.p2p(bs_src, bs_tar, cell_src, cell_tar, rtfmm::LaplaceKernel::KernelType::matrix);
    rtfmm::print_bodies(bs_tar, cell_tar.body.number, cell_tar.body.offset);

    return 0;
}