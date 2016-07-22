#!/bin/bash

export KMP_AFFINITY=compact
export OMP_NUM_THREADS=8

echo "mkl_gemm=["
for size in $(seq 256 256 12000)
do
    ./test_dgemm_mkl.x $size $size $size
    #./test_dgemm_mkl.x $size 100 $size
done
echo "];"
