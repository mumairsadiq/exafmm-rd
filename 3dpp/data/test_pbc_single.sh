#!/bin/bash
#YBATCH -r rtx6000-ada_1
#SBATCH -N 1
#SBATCH -J RTFMM
#SBATCH --time=72:00:00

. /etc/profile.d/modules.sh
nvidia-smi

#./test_pbc -n 24000 -m 128 -P 10 -i 3 --ewald_ksize 21 --use_fft 1
#./test_pbc -P 12 -i 5 --ewald_ksize 51 --use_precompute 1
#./test_pbc -P 10 -i 3 -n 24000 -m 128 --ewald_ksize 21 -v 1
#OPENBLAS_NUM_THREADS=1 ../build/test/test_pbc -i 5 -P 10 -n 24000 -m 128 -v 0
OPENBLAS_NUM_THREADS=1 ../build/test/test_fmm --setting_t 1 -P 20 --th_num 8 -n 24000 -m 128 -v 0
