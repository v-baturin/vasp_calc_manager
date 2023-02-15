#!/bin/bash -l
#SBATCH --account=nn4654k
#SBATCH --job-name=@JOBNAME
#SBATCH --time=0-30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=@NPROCSHARED
#SBATCH --output=slurm.%j.log

# make the program and environment visible to this script
module --quiet purge
module load Gaussian/g16_C.01

# name of input file without extension
input=@INPUT
output=@OUTFILE

# set the heap-size for the job to 20GB
export GAUSS_LFLAGS2="--LindaOptions -s 20000000"
export PGI_FASTMATH_CPU=avx2


# create the temporary folder
export GAUSS_SCRDIR=/cluster/work/users/$USER/$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR

# split large temporary files into smaller parts
lfs setstripe --stripe-count 8 $GAUSS_SCRDIR

# copy input file to temporary folder
cp $SLURM_SUBMIT_DIR/$input $GAUSS_SCRDIR

# run the program
cd $GAUSS_SCRDIR
time g16.ib $input > $output

# copy result files back to submit directory
cp $output $SLURM_SUBMIT_DIR

exit 0
