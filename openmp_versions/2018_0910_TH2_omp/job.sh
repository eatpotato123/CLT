#!/bin/bash
export OMP_NUM_THREADS=8
export OMP_STACKSIZE=1000M
#export KMP_AFFINITY = scatter #none/scatter/compact
yhrun -N 6 -n 16 -c 8 -p nsfc3 a.out > run1.out 
