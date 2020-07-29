#!/bin/sh
#SBATCH -o out
#SBATCH -p lenovo
#SBATCH -J boron
#SBATCH -t 00:00:02
#SBATCH -N 1
#SBATCH -n 8
$GAUSS_EXEDIR/g16 <  B_17_3.gjf  > log