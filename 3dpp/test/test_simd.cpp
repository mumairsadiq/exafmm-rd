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

    real d[4] = {5,6,7,8};
    __m256d d2 = _mm256_broadcast_pd((const __m128d*)(&d[0]));
    print_simd(d2);
    d2 = _mm256_permute_pd(d2,5);
    print_simd(d2);
    __m256d d3 = _mm256_setr_pd(1,2,3,4);
    print_simd(d3);
    __m256d d4 = _mm256_unpacklo_pd(d3,d3);
    print_simd(d4);
    __m256d d5 = _mm256_unpackhi_pd(d3,d3);
    print_simd(d5);
    __m256d d6 = _mm256_addsub_pd(d4,d5);
    print_simd(d6);
    d3 = _mm256_setr_pd(1,2,3,4);
    print_simd(d3);
    d6 = _mm256_permute2f128_pd(d3,d3,1);
    print_simd(d6);
    return 0;
}