#!/bin/bash
#YBATCH -r ad4_1
#SBATCH -N 1
#SBATCH -J RTFMM
#SBATCH --time=24:00:00

. /etc/profile.d/modules.sh
nvidia-smi
for p in $(seq 2 1 20)
do
    OPENBLAS_NUM_THREADS=1 ../build/test/test_fmm --setting_t 1 -P $p --th_num 8 -n 24000 -m 128 -v 0
done
