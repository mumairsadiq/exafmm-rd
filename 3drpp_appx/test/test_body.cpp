#include "type.h"
#include "body.h"

int main(int argc, char* argv[])
{
    std::cout<<"rtfmm_3d_test_body"<<std::endl;

    rtfmm::Bodies3 bs = rtfmm::generate_random_bodies(5, 1);
    rtfmm::print_bodies(bs);

    return 0;
}