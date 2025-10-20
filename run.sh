#!/bin/bash
#SBATCH -N 1                # Total # of nodes
#SBATCH -c 1                # Cores per task requested
#SBATCH -n 32               # Total # of mpi tasks
#SBATCH -t 00:05:00 #(10 min of execution time)
#SBATCH --mem=64GB #(4GB of memory)
#SBATCH -o myjob_%j.o       # Name of stdout output file(%j expands to jobId)
#SBATCH -e myjob_%j.e       # Name of stderr output file(%j expands to jobId)

module load intel impi
srun ./fswm  -t 1 ~/TFM/seqprueba.fna 
echo "done"                 # Write this message on the output file when finnished
 
