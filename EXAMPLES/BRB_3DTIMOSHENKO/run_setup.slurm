#!/bin/bash
#SBATCH --job-name=20X1000
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --time=23:00:00
#SBATCH --mail-user=nb25@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=4
# #SBATCH --exclusive
#SBATCH --partition=commons
#
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST
#
module load MATLAB/2019a
#
matlab -nodisplay -r "fsamp=2^20; ldof=1; DOF='X'; famp=1000; C_transients; quit"
# matlab -nodisplay -r "fsamp=2^20; ldof=8; DOF='Y'; famp=1000; C_transients; quit"
# matlab -nodisplay -r "fsamp=2^20; ldof=6; DOF='Z'; famp=1000; C_transients; quit"

# matlab -nodisplay -r "Ca_savedat('Y'); quit"

