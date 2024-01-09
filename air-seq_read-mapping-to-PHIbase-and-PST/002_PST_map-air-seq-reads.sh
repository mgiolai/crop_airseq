#!/bin/bash
#SBATCH --cpus-per-task=2   # The number of cores
#SBATCH --mem=8G  # The memory per core
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_map.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_map.err
#SBATCH --mail-user=mgiolai.ac@gmail.com
#SBATCH --mail-type=END,FAIL

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/read-data/reads_renamed_merged_fp_cdh_fa

source activate bwa
cd ${references}
#Index PST genome
#srun bwa index Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa


#Map the Air-seq reads to the PST reference
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/bwa_airseq_vs_PST
cd /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/bwa_airseq_vs_PST
for file in $(ls ${reads}/*.fp.cd.fa.gz | grep "airseq_JIC_sample"); do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/.fp.cd.fa.gz//')
	echo ${prefix}
	srun gunzip -c ${file} > ${prefix}.tmp
	source activate bwa
	srun bwa mem -t 4 ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.tmp > ${prefix}.sam
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
	srun rm ${prefix}.tmp ${prefix}.tmp1.sam ${prefix}.tmp2.sam ${prefix}.95pident.95pmatch.bam ${prefix}.95pident.95pmatch.sorted.bam ${prefix}.95pident.95pmatch.sorted.bam.bai
	srun gzip --best ${prefix}.sam
done
