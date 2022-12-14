#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err
#SBATCH --array=0-70

#Don't run this command as an array, needs ~16G ram
#cd /mnt/shared/scratch/mgiolai/air-seq/references
#Index reference
#source activate bwa
#srun bwa index reference.fa
#conda deactivate

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/reads
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_phibase4-12
cd /mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_phibase4-12

ARR=($(ls ${reads}/*.fp.cd.fa.gz))

for file in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/.fp.cd.fa.gz//')
	srun gunzip -c ${file} > ${prefix}.fa
	source activate bwa
	srun bwa mem -t 2 ${references}/reference.fa ${prefix}.fa > ${prefix}.sam
	srun rm ${prefix}.fa
	conda deactivate
	srun filtersam -i 95 -m 95 -o ${prefix}.95pident.95pmatch.sam ${prefix}.sam
	srun gzip --best ${prefix}.sam
	srun gzip --best ${prefix}.95pident.95pmatch.sam
done