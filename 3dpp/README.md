# kernel independent FMM with PBC

## Requirements
+ BLAS
+ LAPACK
+ FFTW3
+ CMAKE

## Compile
```
mkdir build  
cd build  
cmake ..
make  
```

## Executions
+ ./test_fmm : kernel independent FMM without PBC
+ ./test_pbc : kernel independent FMM with PBC

## Parameters
+ -n(int) : number of bodies
+ -c(float) : cycle of periodic images(size of depth0 FMM cell)
+ -P(int) : expansion level
+ -m(int) : cell division threshold
+ -i(int) : depth of periodic images
+ -a(string) : string of algorithms(e.g., 'fe' stands for FMM and ewald-summition, 'fd' stands for FMM and direct-method)
+ -v(bool) : if display detailed information
+ -t(bool) : if timing each step of the computation
+ --use_precompute(bool) : if use precomputation for FMM kernels
+ --seed(int) : random seed for generating bodies
+ --th_num(int) : number of OpenMP threads
+ --dipole_correction(bool) : if use dipole correction for FMM and direct-method

## Typical usags

1. 24000 bodies without PBC, use FMM(P=6), Ewald and direct method
```
./test_fmm -n 24000 -m 128 -P 6 -v 0 -a fde --use_precompute 1
```

2. 24000 bodies with 5-depth peridoic images, use FMM(P=6) and Ewald
```
./test_pbc -n 24000 -m 128 -i 5 -P 10 -v 0 -a fe --use_precompute 1
```