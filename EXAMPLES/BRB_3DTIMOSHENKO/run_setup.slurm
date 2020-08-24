#!/bin/bash
#SBATCH --job-name=IRNZ1000
#SBATCH --ntasks=4
#SBATCH --mem=10G
#SBATCH --time=09:30:00
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
matlab -nodisplay -r "C_transients; quit"

