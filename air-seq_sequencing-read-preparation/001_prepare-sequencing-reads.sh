#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err

#############################
##DEFINE FOLDERS, FILES, OTHER VARIABLES
#############################

reads_raw=/mnt/shared/scratch/mgiolai/air-seq/reads_raw
reads_merged=/mnt/shared/scratch/mgiolai/air-seq/reads
mkdir -p ${reads_merged}

cd ${reads_raw}
for file in $(ls *.fq.gz | xargs -L 1 basename); do
	#Merge and rename lane-split sequencing reads
	filename_short=$(echo ${file}| sed 's/_L001_R1.fastq.gz//')
	filename_output=$(echo ${filename_short} | awk '{split($0,a,"_"); print "airseq_"a[2]"_sample-"a[1]}')
	srun cat ${reads_raw}/${filename_short}_L00*_R1.fastq.gz > ${reads_merged}/${filename_output}.fq.gz
	
	cd ${reads_merged}
	filename_output=$(echo ${file}| sed 's/.fq.gz//')
	##Run fastp - runs at max 16 threads
	source activate fastp
	srun fastp --thread 4 --qualified_quality_phred 20 --length_required 75 --low_complexity_filter --overrepresentation_analysis --adapter_fasta adapters.fasta --html ${filename_output}.html --in1 ${filename_output}.fq.gz --out1 ${filename_output}.fp.fq
	conda deactivate

	#Filter duplicates using cdhit
	source activate cdhit-auxtools
	##Run cd-hit-dup
	srun cd-hit-dup -i ${filename_output}.fp.fq.gz -o ${filename_output}.fp.cd.fq
	source deactivate cdhit-auxtools
	
	#Produce fasta files from cd-hit deduplicated reads
	srun sed -n '1~4s/^@/>/p;2~4p' ${filename_output}.fp.cd.fq > ${filename_output}.fp.cd.fa
done

#Compress all files
srun gzip --best *.fq
srun gzip --best *.fa