#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/reads
ASvsPHI=/mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_phibase4-12

source activate bwa
cd ${references}
#Index PST genome, which has to be downloaded and extracted first
srun bwa index Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa
conda deactivate

#Map the Air-seq reads to the PST reference
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_PST
cd /mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_PST
for file in $(ls ${reads}/*.fp.cd.fa.gz); do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/.fp.cd.fa.gz//')
	echo ${prefix}
	srun gunzip -c ${file} > ${prefix}.tmp
	source activate bwa
	srun bwa mem ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.tmp > ${prefix}.sam
	conda deactivate
	srun rm ${prefix}.tmp
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

#Map the reads that mapped to the PST genome in the PHI-baes database selectively
prefix="phibase.filtered.PST"
echo ${prefix}
srun zcat ${ASvsPHI}/*JIC*95*95*.sam.gz | grep -E "taxid\|27350\|" | grep -v "^@" | awk '{print $1}' > names.txt
srun zcat ${reads}/*JIC*.fa.gz | grep -A1 -f names.txt --no-group-separator > ${prefix}.reads.fa
source activate bwa
srun bwa mem ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.reads.fa > ${prefix}.sam
conda deactivate
srun filtersam -i 95 -m 95 -o ${prefix}.sam ${prefix}.sam
source activate samtools
srun samtools view -S -b ${prefix}.sam > ${prefix}.bam
srun samtools sort ${prefix}.bam -o ${prefix}.sorted.bam
srun samtools index ${prefix}.sorted.bam
srun bcftools mpileup -f ${references}/Puccinia_striiformis.PST-130_1.0.dna.toplevel.fa ${prefix}.sorted.bam > ${prefix}.raw_calls.bcf
srun bcftools call -v -m ${prefix}.raw_calls.bcf > ${prefix}.calls.vcf
srun bcftools filter -i "QUAL>20 && DP>20" ${prefix}.calls.vcf | bcftools view -g ^miss > ${prefix}.filtered.calls.vcf
srun bcftools stats ${prefix}.filtered.calls.vcf | grep "TSTV" > ${prefix}.filtered.calls.TSTV.txt
conda deactivate

#Can we find P. triticina reads?
prefix="phibase.filtered.PT"
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_PT
cd /mnt/shared/scratch/mgiolai/air-seq/bwa_airseq_vs_PT
srun zcat ${ASvsPHI}/*JIC*95*95*.sam.gz | grep -E "taxid\|208348\|" | grep -v "^@" | awk '{print $1}' > names.txt
srun zcat ${reads}/*JIC*.fa.gz | grep -A1 -f names.txt --no-group-separator > ${prefix}.reads.fa