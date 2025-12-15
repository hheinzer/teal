#!/bin/bash
#
# usage: sbatch -n <ntasks> launch.sh <resolution>
#
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --constraint=ib&EPYC_7742
#SBATCH --hint=nomultithread
#SBATCH --exclusive

echo "JobID : $SLURM_JOB_ID"
echo "Nodes : $SLURM_JOB_NODELIST"
echo

export OMP_NUM_THREADS=1
srun --cpu-bind=cores run -c $SLURM_MEM_PER_CPU "$@"
