# kernel independent FMM with PBC
README.md updated : 2024/08/09

## Requirements
+ BLAS
+ LAPACK
+ FFTW3
+ CMAKE

## Compile
```
mkdir build  
cd build  
cmake .. -DCMAKE_BUILD_TYPE=Release
make  
```

## Executions
+ ./test/test_fmm : kernel independent FMM without PBC
+ ./test/test_pbc : kernel independent FMM with PBC
+ ./test/test_reg : kernel independent FMM with PBC and regularization
    - __WARNING__ : currently, enable PBC and Regularization at the same will increase the error

## Parameters
+ -n(int) : number of bodies
+ -c(float) : cycle of periodic images(size of depth0 FMM cell)
+ -P(int) : expansion level
+ -r(float) : size of extension for each regularized cell
+ -m(int) : cell division threshold
+ -i(int) : depth of periodic images
+ -a(string) : string of algorithms(e.g., 'fe' stands for FMM and ewald-summition, 'fd' stands for FMM and direct-method)
+ -v(bool) : if display detailed information
+ -t(bool) : if timing each step of the computation
+ --use_precompute(bool) : if use precomputation for FMM kernels
+ --seed(int) : random seed for generating bodies
+ --th_num(int) : number of OpenMP threads
+ --dipole_correction(bool) : if use dipole correction for FMM and direct-method
+ --setting_t(bool) : set to 1 to quickly match the result of FMM, direct and Ewald

## Typical usages

### 1. 24000 bodies without PBC, use FMM(P=6), Ewald and direct method
```
./test_fmm -n 24000 -m 128 -P 6 -v 0 -a fde
```
result: 
```
-----------------------------FMM vs Direct------------------------------[24000]
L2  (p)  : 5.74209e-08   L2  (f)  : 4.49996e-06   L2  (e)  : 9.02658e-07
Rms (p)  : 2.86722e-05   Rms (f)  : 5.67216e-04
p-energy1 : -2.549454756416e+02
p-energy2 : -2.549452455133e+02
```

### 2. 24000 bodies with 5-depth peridoic images, use FMM(P=10) and Ewald
```
./test_pbc -n 24000 -m 128 -i 5 -P 10 -v 0 -a fe
```
result:  
```
------------------------------FMM vs Ewald------------------------------[24000]
L2  (p)  : 4.48377e+01   L2  (f)  : 2.96419e-07   L2  (e)  : 1.71958e-05
Rms (p)  : 4.99102e+02   Rms (f)  : 3.75112e-05
p-energy1 : 2.773073283084e+01
p-energy2 : 2.773120969004e+01
```

### 3. 1000 bodies without PBC, use FMM(P=6, regularization-size=0.01), and direct method
```
./test_reg -n 1000 -P 6 -m 16 -a fd -r 0.01
```
result:  
```
-----------------------------FMM vs Direct------------------------------[1000]
L2  (p)  : 5.27743e-07   L2  (f)  : 3.19225e-06   L2  (e)  : 9.16410e-07
Rms (p)  : 4.75878e-06   Rms (f)  : 9.54696e-05
p-energy1 : -4.108488165959e+01
p-energy2 : -4.108491931023e+01
```