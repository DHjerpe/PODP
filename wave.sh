#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 00:10
#SBATCH -N 2
#SBATCH -p node

module load gcc openmpi

echo $1 x $1 matrix

echo 1 processor
mpirun -np 1 --map-by core ./wave $1

echo 4 processors
mpirun -np 4 --map-by core ./wave $1

echo 9 processors
mpirun -np 9 --map-by core ./wave $1

echo 12 processors
mpirun -np 12 --map-by core ./wave $1

echo 16 processors
mpirun -np 16 --map-by core ./wave $1

echo 25 processors
mpirun -np 25 --map-by core ./wave $1

echo 36 processors
mpirun -np 36 --map-by core ./wave $1

echo 40 processors
mpirun -np 40 --map-by core ./wave $1
