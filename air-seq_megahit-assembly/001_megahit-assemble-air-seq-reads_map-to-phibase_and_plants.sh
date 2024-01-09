#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --output=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_mhit.out
#SBATCH --error=/mnt/shared/scratch/mgiolai/air-seq/job-outs/slurm-%J_mhit.err
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-70

references=/mnt/shared/scratch/mgiolai/air-seq/references
reads=/mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/read-data/reads_renamed_merged_fp_cdh_fa
mkdir -p /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12

cd ${references}
#Download Chenopodium album from genbank
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/948/465/745/GCA_948465745.1_dcCheAlbu1.1/GCA_948465745.1_dcCheAlbu1.1_genomic.fna.gz
srun gunzip -c GCA_948465745.1_dcCheAlbu1.1_genomic.fna.gz -c | sed "s/^>/>taxid|3559|/g" > chenopodium.fa
#Download Hordeum vulgare from ensembl release-57
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.gz
srun gunzip -c Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.gz -c | sed "s/^>/>taxid|4513|/g" > hordeum.fa
#Download Triticum aestivum from ensembl release-57
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
srun gunzip -c Triticum_aestivum.IWGSC.dna.toplevel.fa.gz -c | sed "s/^>/>taxid|4565|/g" > triticum.fa
#Download Pisum sativum from ensembl release-57
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/pisum_sativum/dna/Pisum_sativum.Pisum_sativum_v1a.dna.toplevel.fa.gz
srun gunzip -c Pisum_sativum.Pisum_sativum_v1a.dna.toplevel.fa.gz -c | sed "s/^>/>taxid|3888|/g" > pisum.fa

srun cat reference.fa chenopodium.fa hordeum.fa triticum.fa pisum.fa > referenceandplants.fa
#Index reference
source activate blast
srun makeblastdb -in referenceandplants.fa -out referenceandplants -dbtype nucl
conda deactivate

ARR=($(ls ${reads}/*.fp.cd.fa.gz))
cd /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12
for file in ${ARR[$SLURM_ARRAY_TASK_ID]}; do
	prefix=$(echo ${file} | xargs -L1 basename | sed 's/.fp.cd.fa.gz//')
	source activate megahit
	srun megahit -t 2 -r ${file} --presets meta-sensitive -o /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12/${prefix}_ms --out-prefix ${prefix}_ms
	conda deactivate
	source activate blast
	#Filter for at least 90% pident, 90% alignment length and an e-value threshold of 1e-50
	srun blastn -task megablast -query /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12/${prefix}_ms/${prefix}_ms.contigs.fa -db ${references}/referenceandplants -max_target_seqs 1 -max_hsps 1 -num_threads 2 -evalue 1e-50 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' | awk '$3 >= 90.0' | awk '$4/$13 >= 0.90' > /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12/${prefix}.blastn.out
	conda deactivate
	srun python /mnt/shared/projects/nhm/mgiolai/air-seq/scripts/20230731_blastn-to-fasta-seqs.py /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12/${prefix}.blastn.out /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12/${prefix}_ms/${prefix}_ms.contigs.fa > /mnt/shared/scratch/mgiolai/air-seq/air-seq_JIC/megahit_airseq_vs_phibase4-12/${prefix}.filtered.fa
done