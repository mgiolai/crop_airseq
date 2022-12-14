#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err
#SBATCH --array=0-70

reads_merged=/mnt/shared/scratch/mgiolai/air-seq/reads
diamondDir=/mnt/shared/scratch/mgiolai/air-seq/diamond
mkdir -p ${diamondDir}

ARR=($(ls ${reads_merged}/*.fp.cd.fq.gz))

cd ${diamondDir}
for reads in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	prefix=$(echo ${reads} | sed 's/.fp.cd.fq.gz//' | xargs -L1 basename)
	echo ${prefix}
	#Align with diamond
	srun /home/mgiola/toolshed/diamond blastx --evalue 1e-10 --mid-sensitive --max-hsps 1 --top 10 --query ${reads} --db /mnt/shared/scratch/mgiolai/air-seq/references/nr_taxid/nr --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --tmpdir /home/cluster/mgiola/scratch --out ${diamondDir}/${prefix}.diamond.out
	#Call taxonomies, calculate alinged reads, extract reads, ignore the taxid 2787854 (i.e. described as 'other entries')
	source activate ete3
	srun python3 /home/mgiola/toolshed/filter-taxonomise-diamond.py ${prefix}.diamond.out.gz ${prefix} 90 0 0 1e-10 0.1 ${reads_merged}/${prefix}.fq.gz ${reads_merged}/${prefix}.fp.fq.gz ${reads_merged}/${prefix}.fp.cd.fq.gz 2787854
	conda deactivate
done

