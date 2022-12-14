#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err
#SBATCH --array=0-55

mkdir -p /mnt/shared/scratch/mgiolai/air-seq/references
cd /mnt/shared/scratch/mgiolai/air-seq/references

#First see in how many files the database is split, this defines the number of jobs for download
ARR=($(seq -w 0 55))
for i in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	echo $i
	wget -c --tries=50 https://ftp.ncbi.nlm.nih.gov/blast/db/nr.${i}.tar.gz
	wget -c --tries=50 https://ftp.ncbi.nlm.nih.gov/blast/db/nr.${i}.tar.gz.md5
	srun tar -xzvf nr.${i}.tar.gz
done

