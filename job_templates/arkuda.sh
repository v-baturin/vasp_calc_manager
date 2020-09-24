#!/bin/sh
#PBS -l nodes=1:ppn=@NPROCSHARED,walltime=10:00:00
#PBS -N @JOBNAME
#PBS -j oe
#PBS -V
cd ${PBS_O_WORKDIR}
module load MPI/impi/2017.0.4 MKL/2017.4.239
$GAUSS_EXEDIR/g16 <  @INPUT  > @OUTFILE
