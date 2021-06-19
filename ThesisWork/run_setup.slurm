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
matlab -nodisplay -r "Nq_pce=$NQPCES; is=dec2base($sis,10)-'0'; is=is(:); pref=sprintf('nlbb_%d',$sis); Nq_pces = ones(1, 7); Nq_pces(is) = Nq_pce; Ir = cell(size(is)); [Ir{:}] = ndgrid(0:Nq_pce-1); Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))'; nxis = Ir*Nq_pce.^((1:length(is))'-1); Irr = zeros(Nq_pce^length(is), 7); Irr(:, is) = Ir; n = $NDES; RQNM_EXPRSURF_PCEFUN(Irr(n, :)+1, nxis(n), Nq_pces, pref); quit"

# matlab -nodisplay -r "Ca_savedat('Y'); quit"



# NEED TO SET VARIABLES EPT, DOF, SMP, WFRC, FMP, REP. Sample:
# sbatch --export=ALL,SMP=18,DOF='X',FMP=0.01,REP=4 --job-name=18X-2 run_setup.slurm

