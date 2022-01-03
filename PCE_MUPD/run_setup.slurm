#!/bin/bash
#SBATCH --job-name=20$-1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=02:30:00
# #SBATCH --mail-user=nb25@rice.edu
# #SBATCH --mail-type=ALL
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=4
# #SBATCH --exclusive
#SBATCH --partition=commons
# #SBATCH --output=/scratch/nb25/octave-jim/ThesisWork/slurmouts/
#
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST
echo "with $@"
#
module load MATLAB/2019a
export MATLAB_LOG_DIR="$SHARED_SCRATCH/nb25/matlab"
matlab -nodisplay -r "mainrun($NDES,$NQUAD); quit;"

# matlab -nodisplay -r "Ca_savedat('Y'); quit"



# NEED TO SET VARIABLES NQUAD and NDES in (1, NQUAD^2)
# sbatch --export=ALL,SMP=18,DOF='X',FMP=0.01,REP=4 --job-name=18X-2 run_setup.slurm

