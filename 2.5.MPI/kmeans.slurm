#!/bin/bash
#SBATCH --job-name=cloudClassifyKMeans
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --mem=MaxMemPerNode

srun python run_k-means_py3.py 1
