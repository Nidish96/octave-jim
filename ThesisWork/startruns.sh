#!/usr/bin/sh

# Coeficient of Friction Study
sis=$1
nqpces=10
ndesigns=$(($nqpces ** ${#sis}))
part=scavenge
if [ $# = 2 ]; 
then
    part=scavenge;
else
    part=$3
fi

for n in `seq 1 $ndesigns`; do
  printf "$n/$ndesigns : "
  sbatch --export=ALL,NDES=$n,NQPCES=$nqpces,sis=$sis --job-name="$2_$n" --partition=$part run_setup.slurm
done

