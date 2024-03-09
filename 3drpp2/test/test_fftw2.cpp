#include <fftw3.h>
#include "type.h"
#include <iostream>

#define N 100

int main()
{
    srand(1234);
    using namespace rtfmm;

    std::vector<real> in(N);
    int nfreq = N / 2 + 1;
    //int nfreq = N;
    int fft_size = nfreq * 2;
    printf("fft_size = %d\n", fft_size);
    std::vector<real> out(fft_size);
    std::vector<real> in2(N);

    for(int i = 0; i < N; i++) 
    {
        in.data()[i] = (double)rand()/RAND_MAX;
    }
    printf("\n Initial Array\n");
    for(int i = 0; i < N; i++) 
    {
        printf("%10f\n", in[i]);
    }

    fftw_plan forward = fftw_plan_dft_r2c_1d(N,in.data(),(fftw_complex*)out.data() ,FFTW_ESTIMATE);
    fftw_execute(forward);

    printf("\n FFT Array\n");
    for(int i = 0; i < nfreq; i++) 
    {
        printf("%10f + %10f i\n", out[i*2+0],out[i*2+1]);
    }

    fftw_plan backward = fftw_plan_dft_c2r_1d(N,(fftw_complex*)out.data(),in2.data() ,FFTW_ESTIMATE);
    fftw_execute(backward);

    real err = 0;
    printf("\n Initial Array\n");
    for(int i = 0; i < N; i++) 
    {
        in2[i] /= N;
        printf("%10f\n", in2[i]);
        err += std::pow(in[i] - in2[i], 2);
    }
    printf("err = %.4f\n", std::sqrt(err / N));


    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
        
    return 0;
}
