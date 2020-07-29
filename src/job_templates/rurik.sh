#!/bin/sh
#SBATCH -o out
#SBATCH -p cpu
#SBATCH -J boron
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
$GAUSS_EXEDIR/g16 <  B_17_3.gjf  > log