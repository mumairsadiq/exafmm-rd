#include "type.h"
#include "fftw3.h"
#include "argument.h"

int main(int argc, char* argv[])
{
    rtfmm::Argument args(argc, argv);
    args.show();

    tbegin(fft);
    typedef rtfmm::real real_t;
    typedef std::vector<real_t> RealVec;

    int N = 2 * args.P - 1;
    int N3 = N * N * N;
    int N_freq = N * N * (N / 2 + 1);
    int n1 = 2 * args.P;

    //int nconv_ = n1 * n1 * n1, nfreq_ = n1 * n1 * (n1 / 2 + 1), NCHILD = 8;
    int nconv_ = N3, nfreq_ = N_freq, NCHILD = 8; n1 = N;
    printf("nconv = %d, nfreq = %d, NCHILD = %d, n1 = %d\n", nconv_, nfreq_, NCHILD, n1);

    RealVec fftw_in(nconv_ * NCHILD);
    RealVec fftw_out(2 * NCHILD * nfreq_);
    int dim[3] = {n1, n1, n1};
    fftw_plan plan = fftw_plan_many_dft_r2c(3, dim, NCHILD,
                                          (real_t*)&fftw_in[0], nullptr, 1, nconv_,
                                          (fftw_complex*)(&fftw_out[0]), nullptr, 1, nfreq_,
                                          FFTW_ESTIMATE);

    #pragma omp parallel for
    for (size_t node_idx=0; node_idx<28; node_idx++) {
      fftw_execute(plan);
    }
    std::cout<<fftw_out[0]<<std::endl;
    tend(fft);
    return 0;
}