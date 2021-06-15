#!/usr/bin/sh

# 1. 1-5067
# 2. 5068-10000
for n in `seq 5068 10000`; do
  sbatch --export=ALL,NDES=$n --job-name="D_$n" run_setup.slurm
done
