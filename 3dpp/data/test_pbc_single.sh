#!/bin/bash
#YBATCH -r 6000_1
#SBATCH -N 1
#SBATCH -J RTFMM
#SBATCH --time=72:00:00

. /etc/profile.d/modules.sh
nvidia-smi

#./test_pbc -n 24000 -m 128 -P 10 -i 3 --ewald_ksize 21 --use_fft 1
#./test_pbc -P 12 -i 5 --ewald_ksize 51 --use_precompute 1
./test_pbc -P 10 -i 3 -n 24000 -m 128 --ewald_ksize 21 -v 1
