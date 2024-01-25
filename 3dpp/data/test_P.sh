#!/bin/bash
#YBATCH -r rtx6000-ada_1
#SBATCH -N 1
#SBATCH -J RTFMM
#SBATCH --time=24:00:00

. /etc/profile.d/modules.sh
nvidia-smi
for p in $(seq 2 1 20)
do
    OPENBLAS_NUM_THREADS=1 ./test_fmm --zero_netcharge 0 --dipole_correction 0 --divide_4pi 1 -P $p --th_num 8 -n 24000 -m 128 --use_simd 1 -v 0
done
