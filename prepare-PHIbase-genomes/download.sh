#!/bin/bash

##################################
#Make dirs and change to folder
##################################
mkdir download
cp download.txt download
cd download
cat download.txt | awk -v FS='\t' '{print $5}' | sed $'s/\r$//' | xargs -n 1 -P 10 wget -c 
if [ -f md5checksums.all.txt ]; then rm md5checksums.all.txt; fi
cat download.txt | awk -v FS='\t' '{print $4}' | sed $'s/\r$//' | xargs -n 1 -P 10 wget -c  -O - >> md5checksums.all.txt

#Files can be corrupted from download; recover here and then with a MANUAL FINAL CHECK
files=$(cat download.txt | awk -v FS='\t' '{print $5}' | sed $'s/\r$//' | xargs -L 1 basename)
for file in ${files[@]}; do
	if [ -f ${file} ]; then
		md5=$(md5sum ${file} | awk '{print $1}')
		if grep -q ${md5} md5checksums.all.txt; then
			echo ${file}": OK"
		else
			echo ${file}": NOT OK"
			grep ${file} download.txt | awk -v FS='\t' '{print $5}' | sed $'s/\r$//' | xargs -n 1 wget -c
		fi
	fi
done

#Check taxids before modifying the fasta headers
if [ -f taxids.txt ]; then rm taxids.txt; fi
for file in ${files[@]}; do
	taxid=$(grep ${file} download.txt | awk -v FS='\t' '{print $2}')
	echo ${taxid} >> taxids.txt
done

#Modify fasta headers by assigning the taxid number in front of the description and cat to single fasta file
if [ -f reference.fa ]; then rm reference.fa; fi
for file in ${files[@]}; do
	echo ${file}
	taxid=$(grep ${file} download.txt | awk -v FS='\t' '{print $2}')
	gunzip -c ${file} | sed "s/^>/>taxid|${taxid}|/g" >> reference.fa
done

#Download and add human T2T genome to reference
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gunzip -c GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz | sed "s/^>/>taxid|9606|/g" >> reference.fa

#Count genome sizes by taxid
if [ -f reference_genomesizes.txt ]; then rm reference_genomesizes.txt; fi
for file in ${files[@]}; do echo ${file} && taxid=$(grep ${file} download.txt | awk -v FS='\t' '{print $2}') && nt=$(gunzip -c ${file} | grep -v "$>" | wc -c) && echo -e ${taxid}'\t'${nt}'\t'${file} >> reference_genomesizes.txt; done

#Add human reference genome size
nt=$(gunzip -c GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz | grep -v "$>" | wc -c)
echo -e '9606\t'${nt}'\tGCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz' >> reference_genomesizes.txt
