#!/bin/bash

#Submit this script with: sbatch thefilename
# note: change the memory, threads, wall, etc

#SBATCH -t 24:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 272   # total number of processor cores in this job; each node has 272 cores
#SBATCH -J "Dconcat"   # job name
#SBATCH --mem=350 # how much memory you need; each box has ~350G
#SBATCH --output=slurm-Dconcat.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load raxml-knl/8.2.10
raxmlHPC-PTHREADS-AVX2 -f a -k -m GTRGAMMA -o F1_1 -p 25632 -T 272 -s all.phy -x 451842 -n all.GTRGAMMA.raxml -N 10000
