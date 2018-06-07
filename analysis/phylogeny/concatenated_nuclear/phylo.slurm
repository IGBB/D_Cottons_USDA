#!/bin/bash

#Submit this script with: sbatch thefilename
# note: change the memory, threads, wall, etc

#SBATCH -t 148:00:00   # walltime
#SBATCH -N 3   # number of nodes in this job
#SBATCH -J "raxml"   # job name
#SBATCH --mem=340G # how much memory you need; each box has ~340G
#SBATCH --output=raxmlungapped.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

raxmlHPC-MPI-SSE3 -f a -k -m GTRGAMMA -o F1_1 -p 25632 -s small.ungapped.fasta -x 451842 -n ungapped.GTRGAMMA.raxml -N 1000000
