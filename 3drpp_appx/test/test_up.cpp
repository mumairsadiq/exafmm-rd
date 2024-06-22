#include "type.h"
#include "align.h"
#include <iostream>
#include <memory>
#include <fftw3.h>
#include <cstring>
#include "argument.h"
#include <omp.h>

const int MEM_ALIGN = 64;
const int CACHE_SIZE = 512;
const int NCHILD = 8;

typedef double real_t;
typedef rtfmm::AlignedAllocator<real_t, MEM_ALIGN> AlignAllocator;
typedef std::vector<real_t> RealVec;
typedef std::vector<real_t, AlignAllocator> AlignedVec;

int nsurf = 2168;
int nconv = 64000;
int nfreq = 33600;

RealVec surface(int p, real_t r0, int level, real_t * c, real_t alpha) {
    int n = 6*(p-1)*(p-1) + 2;
    RealVec coord(n*3);
    coord[0] = -1.0;
    coord[1] = -1.0;
    coord[2] = -1.0;
    int count = 1;
    for (int i=0; i<p-1; i++) {
      for (int j=0; j<p-1; j++) {
        coord[count*3  ] = -1.0;
        coord[count*3+1] = (2.0*(i+1)-p+1) / (p-1);
        coord[count*3+2] = (2.0*j-p+1) / (p-1);
        count++;
      }
    }
    for (int i=0; i<p-1; i++) {
      for (int j=0; j<p-1; j++) {
        coord[count*3  ] = (2.0*i-p+1) / (p-1);
        coord[count*3+1] = -1.0;
        coord[count*3+2] = (2.0*(j+1)-p+1) / (p-1);
        count++;
      }
    }
    for (int i=0; i<p-1; i++) {
      for (int j=0; j<p-1; j++) {
        coord[count*3  ] = (2.0*(i+1)-p+1) / (p-1);
        coord[count*3+1] = (2.0*j-p+1) / (p-1);
        coord[count*3+2] = -1.0;
        count++;
      }
    }
    for (int i=0; i<(n/2)*3; i++) {
      coord[count*3+i] = -coord[i];
    }
    real_t r = r0 * powf(0.5, level);
    real_t b = alpha * r;
    for (int i=0; i<n; i++){
      coord[i*3+0] = coord[i*3+0]*b + c[0];
      coord[i*3+1] = coord[i*3+1]*b + c[1];
      coord[i*3+2] = coord[i*3+2]*b + c[2];
    }
    return coord;
  }

std::vector<int> generate_surf2conv_up(int p) {
    int n1 = 2*p;
    real_t c[3];
    for (int d=0; d<3; d++) c[d] = 0.5*(p-1);
    RealVec surf = surface(p, 0.5, 0, c, real_t(p-1));
    std::vector<int> map(6*(p-1)*(p-1)+2);
    for (size_t i=0; i<map.size(); i++) {
      map[i] = (int)(p-1-surf[i*3])
             + ((int)(p-1-surf[i*3+1])) * n1
             + ((int)(p-1-surf[i*3+2])) * n1 * n1;
    }
    return map;
  }

void fft_up_equiv(RealVec& all_up_equiv, AlignedVec& fft_in, int p) {
    TIME_BEGIN(up_prepare);
    int& nsurf_ = nsurf;
    int& nconv_ = nconv;
    int& nfreq_ = nfreq;
    
    int n1 = 2 * p;
    auto map = generate_surf2conv_up(p);

    size_t fft_size = 2 * NCHILD * nfreq_;
    printf("nsurf_ = %d, nconv_ = %d, nfreq_ = %d, fft_size = %ld\n",nsurf_,nconv_,nfreq_,fft_size);
    AlignedVec fftw_in(nconv_ * NCHILD);
    AlignedVec fftw_out(fft_size);
    //RealVec fftw_in(nconv_ * NCHILD);
    //RealVec fftw_out(fft_size);
    int dim[3] = {n1, n1, n1};
    fftw_plan plan = fftw_plan_many_dft_r2c(3, dim, NCHILD,
                                          (real_t*)&fftw_in[0], nullptr, 1, nconv_,
                                          (fftw_complex*)(&fftw_out[0]), nullptr, 1, nfreq_,
                                          FFTW_ESTIMATE);
    TIME_END(up_prepare);
    TIME_BEGIN(up_core);
    #pragma omp parallel for
    for (size_t node_idx=0; node_idx<72; node_idx++) {
      RealVec buffer(fft_size, 0);
      //real_t* up_equiv = &all_up_equiv[0];  // offset ptr of node's 8 child's upward_equiv in all_up_equiv, size=8*nsurf_
      real_t* up_equiv_f = &fft_in[fft_size*node_idx]; // offset ptr of node_idx in fft_in vector, size=fft_size
      //std::memset(up_equiv_f, 0, fft_size*sizeof(real_t));  // initialize fft_in to 0
      //for(int i = 0; i < fft_size; i++) up_equiv_f[i] = 0;
      /*for (int k=0; k<nsurf_; k++) {
        size_t idx = map[k];
        for (int j=0; j<NCHILD; j++)
          up_equiv_f[idx+j*nconv_] = 0;
      }*/
      fftw_execute_dft_r2c(plan, up_equiv_f, (fftw_complex*)&buffer[0]);
      /*for (int k=0; k<nfreq_; k++) {
        for (int j=0; j<NCHILD; j++) {
          up_equiv_f[2*(NCHILD*k+j)+0] = buffer[2*(nfreq_*j+k)+0];
          up_equiv_f[2*(NCHILD*k+j)+1] = buffer[2*(nfreq_*j+k)+1];
        }
      }*/
      //fftw_execute(plan);
    }
    TIME_END(up_core);
    fftw_destroy_plan(plan);
  }

int main(int argc, char* argv[])
{
    rtfmm::title("test_up");
    rtfmm::Argument args(argc, argv);
    args.show();

    omp_set_num_threads(args.th_num);
    if(rtfmm::verbose) printf("# of threads = %d\n", omp_get_max_threads());

    int nnodes = 585;
    int nsurf_ = 2168;
    size_t fft_size = 2 * NCHILD * nfreq;
    printf("fft_size = %ld\n", fft_size);
    // allocate memory
    std::vector<real_t> all_up_equiv, all_dn_equiv;
    all_up_equiv.reserve(nnodes*nsurf_);   // use reserve() to avoid the overhead of calling constructor
    all_dn_equiv.reserve(nnodes*nsurf_);   // use pointer instead of iterator to access elements 
    AlignedVec fft_in, fft_out;
    fft_in.reserve(72*fft_size);
    fft_out.reserve(72*fft_size);
      /*for (int i=0; i<nnodes; i++) {
        for (int j=0; j<nsurf_; j++) {
          all_up_equiv[i*nsurf_+j] = i;
          all_dn_equiv[i*nsurf_+j] = j;
        }
      }*/

    TIME_BEGIN(fft_up);
    fft_up_equiv(all_up_equiv, fft_in, args.P);
    TIME_END(fft_up);

    std::cout<<fft_in[0]<<std::endl;
    return 0;
}