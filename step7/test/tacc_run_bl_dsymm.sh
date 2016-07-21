#!/bin/bash
#SBATCH -J bl_dgemm_job
#SBATCH -o bl_dgemm_output-%j.txt
#SBATCH -p gpu
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A CompEdu
export OMP_NUM_THREADS=8
export BLISGEMM_IC_NT=8
export KMP_AFFINITY=compact,verbose

ibrun tacc_affinity run_bl_dsymm.sh
