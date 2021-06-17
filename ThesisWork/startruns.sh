#!/usr/bin/sh

# Coeficient of Friction Study
sis=12346
nqpces=7

for n in `seq 1 100`; do
  printf "$n : "
  sbatch --export=ALL,NDES=$n,NQPCES=$nqpces,sis=$sis --job-name="ROT_$n" --partition=scavenge run_setup.slurm
done

