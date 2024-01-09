#!/bin/bash
#SBATCH --cpus-per-task=6   # The number of cores
#SBATCH --mem=24G  # The memory per core
#SBATCH --partition=medium
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_bwamap.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_bwamap.err
#SBATCH --mail-user=mgiolai.ac@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-70

#2cores 16G
cd /mnt/shared/scratch/mgiolai/air-seq/references/
##Index reference
#source activate bwa
#srun bwa index reference.fa
#conda deactivate

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/read-data/reads_renamed_merged_fp_cdh_fa
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/bwa_airseq_vs_phibase4-12
cd /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/bwa_airseq_vs_phibase4-12

ARR=($(ls ${reads}/*.fp.cd.fa.gz))

for file in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/.fp.cd.fa.gz//')
	srun gunzip -c ${file} > ${prefix}.fa
	source activate bwa
	srun bwa mem -t 12 ${references}/reference.fa ${prefix}.fa > ${prefix}.sam
	srun rm ${prefix}.fa
	conda deactivate
	srun filtersam -i 95 -o ${prefix}.tmp1.sam ${prefix}.sam
	srun filtersam -m 95 -o ${prefix}.tmp2.sam ${prefix}.tmp1.sam
	source activate bwa
	srun cat ${prefix}.tmp2.sam | samtools view -q 1 > ${prefix}.95pident.95pmatch.sam 
	conda deactivate
	srun rm ${prefix}.tmp1.sam
	srun rm ${prefix}.tmp2.sam
	srun gzip --best ${prefix}.sam
	srun gzip --best ${prefix}.95pident.95pmatch.sam
done
