#!/bin/sh
#PBS -l nodes=1:ppn=10,walltime=0:01:00
#PBS -N  Cd14Se2i1_
#PBS -j oe
#PBS -V
cd ${PBS_O_WORKDIR}
module load MPI/impi/2017.0.4 MKL/2017.4.239
$GAUSS_EXEDIR/g16 <  Cd14Se2.gjf  > log