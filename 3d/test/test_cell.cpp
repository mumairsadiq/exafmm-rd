#include "type.h"
#include "body.h"
#include "cell.h"

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3d_test_cell"<<std::endl;

    rtfmm::real r = 1.0;
    int num_body = 5;

    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(num_body, r);

    rtfmm::Cell3 cell;
    cell.idx = 0;
    cell.depth = 0;
    cell.r = r;
    cell.x = {0,0,0};
    cell.child = {0,0};
    cell.body = {0,num_body};

    rtfmm::print_bodies(bs, cell.body.number, cell.body.offset);

    return 0;
}