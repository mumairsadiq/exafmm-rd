#!/bin/bash
#YBATCH -r rtx6000-ada_1
#SBATCH -N 1
#SBATCH -J RTFMM
#SBATCH --time=24:00:00

. /etc/profile.d/modules.sh
nvidia-smi
lscpu
ip a
echo machine-id=$(cat /etc/machine-id)

for p in $(seq 4 1 8)
do
    for i in $(seq 0 1 9)
    do
        OPENBLAS_NUM_THREADS=1 ../build/test/test_pbc -P $p --th_num 8 -n 24000 -m 128 -v 0 -i $i -a fe
    done
done
