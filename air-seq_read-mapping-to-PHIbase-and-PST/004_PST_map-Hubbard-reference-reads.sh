#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err
#SBATCH --array=0-66

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/reads/pst_reference

ARR=($(ls ${reads}/*_1.fastq.gz))

cd air-seq/air-seq_JIC/bwa_airseq_vs_PST
for file in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/_1.fastq.gz//' )
	echo ${prefix}
	source activate fastp
	srun fastp -i ${reads}/${prefix}_1.fastq.gz -I ${reads}/${prefix}_2.fastq.gz -o ${prefix}.R1.fq -O ${prefix}.R2.fq
	conda deactivate
	source activate bwa
	srun bwa mem ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.R1.fq ${prefix}.R2.fq > ${prefix}.sam
	conda deactivate
	srun rm ${prefix}.R1.fq ${prefix}.R2.fq
	srun filtersam -i 95 -m 95 -o ${prefix}.95pident.95pmatch.sam ${prefix}.sam
	source activate samtools
	srun samtools view -S -b ${prefix}.95pident.95pmatch.sam > ${prefix}.95pident.95pmatch.bam
	srun samtools sort ${prefix}.95pident.95pmatch.bam -o ${prefix}.95pident.95pmatch.sorted.bam
	srun samtools index ${prefix}.95pident.95pmatch.sorted.bam
	srun bcftools mpileup -f ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.95pident.95pmatch.sorted.bam > ${prefix}.95pident.95pmatch.raw_calls.bcf
	srun bcftools call -v -m ${prefix}.95pident.95pmatch.raw_calls.bcf > ${prefix}.95pident.95pmatch.calls.vcf
	srun bcftools filter -i "QUAL>20 && DP>20" ${prefix}.95pident.95pmatch.calls.vcf | bcftools view -g ^miss > ${prefix}.95pident.95pmatch.filtered.calls.vcf
	srun bcftools stats ${prefix}.95pident.95pmatch.filtered.calls.vcf | grep "TSTV" > ${prefix}.95pident.95pmatch.filtered.calls.TSTV.txt
	conda deactivate
done