#!/bin/bash
#SBATCH --cpus-per-task=2   # The number of cores
#SBATCH --mem=8G  # The memory per core
#SBATCH --partition=medium
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_map.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_map.err
#SBATCH --mail-user=mgiolai.ac@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-66

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/read-data/pst_reference

ARR=($(ls ${reads}/*_1.fastq.gz))

mkdir -p /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/bwa_airseq_vs_PST
cd /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/bwa_airseq_vs_PST
for file in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/_1.fastq.gz//' )
	echo ${prefix}
	source activate fastp
	srun fastp --thread 4 --qualified_quality_phred 20 --length_required 75 --low_complexity_filter -i ${reads}/${prefix}_1.fastq.gz -I ${reads}/${prefix}_2.fastq.gz -o ${prefix}.R1.fq -O ${prefix}.R2.fq -h ${prefix}.report.html
	conda deactivate
	source activate bwa
	srun bwa mem -t 4 ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.R1.fq ${prefix}.R2.fq > ${prefix}.sam
	conda deactivate
	srun filtersam -i 95 -o ${prefix}.tmp1.sam ${prefix}.sam
	srun filtersam -m 95 -o ${prefix}.tmp2.sam ${prefix}.tmp1.sam
	source activate bwa
	srun samtools view -b -q 1 ${prefix}.tmp2.sam > ${prefix}.95pident.95pmatch.bam
	srun samtools sort ${prefix}.95pident.95pmatch.bam -o ${prefix}.95pident.95pmatch.sorted.bam
	srun samtools index ${prefix}.95pident.95pmatch.sorted.bam
	srun bcftools mpileup -f ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.95pident.95pmatch.sorted.bam > ${prefix}.95pident.95pmatch.raw_calls.bcf
	srun bcftools call -v -m ${prefix}.95pident.95pmatch.raw_calls.bcf > ${prefix}.95pident.95pmatch.calls.vcf
	srun bcftools filter -i "QUAL>20 && DP>20" ${prefix}.95pident.95pmatch.calls.vcf | bcftools view -g ^miss > ${prefix}.95pident.95pmatch.filtered.calls.vcf
	srun bcftools stats ${prefix}.95pident.95pmatch.filtered.calls.vcf | grep "TSTV" > ${prefix}.95pident.95pmatch.filtered.calls.TSTV.txt
	conda deactivate
	srun rm ${prefix}.R1.fq ${prefix}.R2.fq ${prefix}.tmp1.sam ${prefix}.tmp2.sam ${prefix}.95pident.95pmatch.bam ${prefix}.95pident.95pmatch.sorted.bam ${prefix}.95pident.95pmatch.sorted.bam.bai
	srun pigz --best -p4 ${prefix}.sam
done
