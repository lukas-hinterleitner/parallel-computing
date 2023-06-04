#!/bin/bash
#SBATCH -N 6
#SBATCH --ntasks 96
#SBATCH -t 4000
mpirun ./a2-mpi --m 1152 --n 1152 --epsilon 0.01 --max-iterations 1000 --no-verify