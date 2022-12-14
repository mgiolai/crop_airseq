#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=long
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J.err

cd /mnt/shared/scratch/mgiolai/air-seq/references
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/references/nr_taxid

#Add the taxid to the header of each fasta sequence in the nr database, separate the taxid and the description with &&& for later processing
srun /home/mgiola/toolshed/ncbi-blast-2.12.0+/bin/blastdbcmd -entry all -db nr -outfmt ">%T&&&%a &&& %s" | sed 's/ &&& /\n/g' >> /mnt/shared/scratch/mgiolai/air-seq/references/nr_taxid/nr.taxid.amended.fa
#Make DIAMOND database
cd /mnt/shared/scratch/mgiolai/air-seq/references/nr_taxid
srun /home/mgiola/toolshed/toolshed/diamond makedb --in nr.taxid.amended.fa -d nr