#include "type.h"

int main(int argc, char* argv[])
{
    using namespace rtfmm;
    title("rtfmm_3dpp_test_SIMD");
    vec4d v1(1,2,3,4);
    vec4d v2(1,2,3,5);
    std::cout<<v1<<std::endl;
    std::cout<<v2<<std::endl;
    std::cout<<(v1==v2)<<std::endl;
    std::cout<<(v1!=v2)<<std::endl;
    std::cout<<(v1+v2)<<std::endl;
    std::cout<<(v1-v2)<<std::endl;
    std::cout<<(v1*v2)<<std::endl;
    std::cout<<(v1/v2)<<std::endl;
    std::cout<<(-v1)<<std::endl;
    std::cout<<v1.norm()<<std::endl;
    std::cout<<v1.r()<<std::endl;
    std::cout<<v1.sum()<<std::endl;
    v1 += 1;
    std::cout<<v1<<std::endl;
    v1 -= 1;
    std::cout<<v1<<std::endl;
    v1 *= 2;
    std::cout<<v1<<std::endl;
    v1 /= 4;
    std::cout<<v1<<std::endl;
    std::cout<<v1.sqrt()<<std::endl;

    vec4d v3 = v1;
    v3 &= v2 <= 2;
    std::cout<<v3<<std::endl;
    return 0;
}