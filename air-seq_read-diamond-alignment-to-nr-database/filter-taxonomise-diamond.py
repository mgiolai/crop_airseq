#!/usr/bin/env python

import os
import sys
import gzip
import time
from ete3 import NCBITaxa
from Bio import SeqIO
from collections import Counter
import progressbar
start_time = time.time()

#DIAMOND output columns that need to be present (python numbered)
#0 qseqid
#1 sseqid
#2 pident
#3 length
#4 mismatch
#5 gapopen
#6 qstart
#7 qend
#8 sstart
#9 send
#10 evalue
#11 bitscore
#12 qlen
#13 slen

#############################
#Input variables
#############################

#INPUT variabels that need to be present
diamondFile=sys.argv[1]
filePrefix=str(sys.argv[2])
pBSThreshold=float(sys.argv[3])
pnmatchThreshold=float(sys.argv[4])
pidentThreshold=float(sys.argv[5])
evalueThreshold=float(sys.argv[6])
LCAreadSupportThreshold=float(sys.argv[7])
rawReadFastqFile=sys.argv[8]
qualityTrimFastqFile=sys.argv[9]
DedupFastqFile=sys.argv[10]
ignoreIDs=set(str(sys.argv[11]).split(','))

print("File: "+str(diamondFile),flush=True)
print("Prefix: "+str(filePrefix),flush=True)
print("Bitscore cutoff: "+str(pBSThreshold),flush=True)
print("Percent nmatch cutoff: "+str(pnmatchThreshold),flush=True)
print("Percent identity cutoff: "+str(pidentThreshold),flush=True)
print("e-value threshold: "+str(evalueThreshold),flush=True)
print("Minimal percent read support / taxon: "+str(LCAreadSupportThreshold),flush=True)
print("Raw read fastq file: "+str(rawReadFastqFile),flush=True)
print("Quality trimmed read fastq file: "+str(qualityTrimFastqFile),flush=True)
print("Aligned read fastq file: "+str(DedupFastqFile),flush=True)

#############################
#Output files and filenames
#############################

filenameOut1=str(filePrefix)+'_reads_nonpass.fa'
filenameOut2=str(filePrefix)+'_readnumbers.txt'

filenameOut3='001_1_{}_{}BS-{}aln_{}idt_{}ev_NArm_SK_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut4='001_2_{}_{}BS-{}aln_{}idt_{}ev_NArm_PHY_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut5='001_3_{}_{}BS-{}aln_{}idt_{}ev_NArm_CL_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut6='001_4_{}_{}BS-{}aln_{}idt_{}ev_NArm_OR_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut7='001_5_{}_{}BS-{}aln_{}idt_{}ev_NArm_FA_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut8='001_6_{}_{}BS-{}aln_{}idt_{}ev_NArm_GE_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut9='001_7_{}_{}BS-{}aln_{}idt_{}ev_NArm_SP_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))
filenameOut10='001_8_{}_{}BS-{}aln_{}idt_{}ev_NArm_ST_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold))

filenameOut11='002_1_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-SK_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut12='002_2_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-PHY_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut13='002_3_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-CL_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut14='002_4_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-OR_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut15='002_5_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-FA_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut16='002_6_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-GE_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut17='002_7_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-SP_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))
filenameOut18='002_8_{}_{}BS-{}aln_{}idt_{}ev_NArm_{}pMin-ST_cts.txt'.format(filePrefix, str(pBSThreshold), str(pnmatchThreshold), str(pidentThreshold), str(evalueThreshold), str(LCAreadSupportThreshold))

#############################
#Update NCBI taxonomy database
#############################
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

#############################
#Definitions
#############################

def get_desired_ranks(taxid, desired_ranks):
	#Taxonomise based on a taxid using the ete3 toolkit default commands
	taxonomy=[]
	try:
		lineage = ncbi.get_lineage(int(taxid))
		#names = ncbi.get_taxid_translator(lineage)
		lineage2ranks = ncbi.get_rank(ncbi.get_taxid_translator(lineage))
		ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
		ranks={'{}_id'.format(rank): ranks2lineage.get(rank, 'NA') for rank in desired_ranks}
		if ranks != None:
			for key, rank in ranks.items():
				if rank != 'NA':
					taxonomy.append(list(ncbi.get_taxid_translator([rank]).values())[0])
				else:
					taxonomy.append('NA')
		else:
			taxonomy.append('NONE')
		return(taxonomy)
	except:
		taxonomy.append('NONE')
		return(taxonomy)

def countreads(gzipFastqFile):
	ct=sum(1 for seq in SeqIO.parse(gzip.open(gzipFastqFile, 'rt'), 'fastq'))
	return(ct)

def addQnameToTaxDict(qname,taxonomy,dictQnameToTaxonomy):
	#Prepare a taxonomy dictionary for each read with the specified ranks from superkingdom to strain
	if qname not in dictQnameToTaxonomy.keys():
		#Create dictQnameToTaxonomy
		dictQnameToTaxonomy[qname]={}
		#Add nested sets
		dictQnameToTaxonomy[qname]['superkingdom']=set()
		dictQnameToTaxonomy[qname]['phylum']=set()
		dictQnameToTaxonomy[qname]['class']=set()
		dictQnameToTaxonomy[qname]['order']=set()
		dictQnameToTaxonomy[qname]['family']=set()
		dictQnameToTaxonomy[qname]['genus']=set()
		dictQnameToTaxonomy[qname]['species']=set()
		dictQnameToTaxonomy[qname]['strain']=set()
		#Add taxonomy to nested sets
		dictQnameToTaxonomy[qname]['superkingdom'].add(taxonomy[0])
		dictQnameToTaxonomy[qname]['phylum'].add(taxonomy[1])
		dictQnameToTaxonomy[qname]['class'].add(taxonomy[2])
		dictQnameToTaxonomy[qname]['order'].add(taxonomy[3])
		dictQnameToTaxonomy[qname]['family'].add(taxonomy[4])
		dictQnameToTaxonomy[qname]['genus'].add(taxonomy[5])
		dictQnameToTaxonomy[qname]['species'].add(taxonomy[6])
		dictQnameToTaxonomy[qname]['strain'].add(taxonomy[7])
	else:
		#Add taxonomy to nested sets
		dictQnameToTaxonomy[qname]['superkingdom'].add(taxonomy[0])
		dictQnameToTaxonomy[qname]['phylum'].add(taxonomy[1])
		dictQnameToTaxonomy[qname]['class'].add(taxonomy[2])
		dictQnameToTaxonomy[qname]['order'].add(taxonomy[3])
		dictQnameToTaxonomy[qname]['family'].add(taxonomy[4])
		dictQnameToTaxonomy[qname]['genus'].add(taxonomy[5])
		dictQnameToTaxonomy[qname]['species'].add(taxonomy[6])
		dictQnameToTaxonomy[qname]['strain'].add(taxonomy[7])
	return(dictQnameToTaxonomy)

def findTaxonomies(dictQnameToTaxonomy,taxonomyFilter):
	#Find the taxonomies which are uniquely assigned to a read:
	#1. Specify taxonomic levels of interest as a list using the taxonomyFilter variable
	#2. For each read count if the read hits the specified taxonomic levels uniquely (i.e. not multiple superkingdoms, classes, orders, species, etc. just one)
	#3. If a read hits the specified taxonomic levels uniquely add the read and the taxonomic levels to the filteredDictionary
	#4. Last take all found taxonomies in the filteredDictionary and load them to the uniqueTaxonomyfilteredDictionary.
	# The taxonomies are only added once to the uniqueTaxonomyfilteredDictionary and later on utilised to combine and output the described ranks with their corresponding counts
	filteredDictionary={}
	for qname, taxonomy in dictQnameToTaxonomy.items():
		ct=0
		for rank in taxonomyFilter:
			if len(taxonomy[rank])==1:
				ct+=1
		if ct == len(taxonomyFilter):
			for rank, name in taxonomy.items():
				if rank in taxonomyFilter:
					name=list(name)[0]
					if qname not in filteredDictionary.keys():
						filteredDictionary[qname]={}
						filteredDictionary[qname][rank]=set()
						filteredDictionary[qname][rank].add(name)
					else:
						filteredDictionary[qname][rank]=set()
						filteredDictionary[qname][rank].add(name)
	#Subset so that only unique taxonomies are present in the dictionary to avoid doubled taxonomy printing in the output
	uniqueTaxonomyfilteredDictionary={}
	for qname, taxonomy in filteredDictionary.items():
		if taxonomy not in uniqueTaxonomyfilteredDictionary.values():
			uniqueTaxonomyfilteredDictionary[qname]=taxonomy
	return(uniqueTaxonomyfilteredDictionary)


#############################
#Filtering and output
#############################

#Define the NCBI taxids and descendant IDs that will be ignored
print('Obtaining taxids and their descendants that will be ignored by user definition',flush=True)
taxidsIgnore=set()
for ID in ignoreIDs:
	taxidsIgnore.add(ID)
	for descendantID in ncbi.get_descendant_taxa(ID):
		taxidsIgnore.add(descendantID)

#Determine read numbers
print('Counting number of raw reads',flush=True)
nreadsAll=countreads(rawReadFastqFile)

print('Counting number of quality trimmed reads',flush=True)
nreadsTrim=countreads(qualityTrimFastqFile)

print('Counting number of deduplicated reads',flush=True)
nreadsDedup=countreads(DedupFastqFile)

#Determine highest bitscore per read
print('Determining highest bitscore per read',flush=True)
#Calculate total number of lines in diamond file for setting up progress bar
n=sum(1 for i in gzip.open(diamondFile, 'rt'))
mark=round(n/10)
ct=0
#Extract bitscores
dictBS={}
with gzip.open(diamondFile, 'rt') as diamondOut:
	with progressbar.ProgressBar(max_value=n) as pbar:
		pbar.update(ct)
		for line in diamondOut:
			line=line.replace('\n','').split('\t')
			qname, score=str(line[0]), float(line[11])
			if qname not in dictBS.keys() or score > dictBS[qname]:
				dictBS[qname]=score
			ct+=1
			if ct%mark==0:
				pbar.update(ct)
diamondOut.close()

##Calculate bitscore %-cutoff
print('Calculating bitscore cutoffs',flush=True)
f=float(pBSThreshold)/100
dictBS.update({qname:f*dictBS[qname] for qname in dictBS.keys()})

#Filter for %-bitscore per read, %-identity, %-aligned bases cutoff, evalue and taxonomise reads
print('Filtering bitscore, percent-match, percent-identity, evalue, taxonomising',flush=True)
setReadsAlignAfterFiltering=set()
setReadsAlignAfterFilteringTaxonomised=set()
dictQnameToTaxonomy={}
dictTaxidToTaxonomy={}
ct=0
with gzip.open(diamondFile, 'rt') as diamondOut:
	with progressbar.ProgressBar(max_value=n) as pbar:
		pbar.update(ct)
		for line in diamondOut:
			line=line.replace('\n','').split('\t')
			qname, score, pnmatch, pident, evalue=str(line[0]), float(line[11]), 100*(float(line[12])-float(line[4]))/float(line[12]), float(line[2]), float(line[10])
			if score >= float(dictBS[qname]) and pnmatch >= pnmatchThreshold and pident >= pidentThreshold and evalue <= evalueThreshold:
				#Add to read set that passed filtering
				setReadsAlignAfterFiltering.add(qname)
				#Extract taxid from subject name
				taxid=int(str(line[1]).split('&&&')[0])
				#IDs from the taxid 'other entries' are not processed
				if taxid not in taxidsIgnore:
					#Obtain taxonomy, if already known pull from dictTaxidToTaxonomy otherwise predict and add to dictTaxidToTaxonomy
					if taxid in dictTaxidToTaxonomy.keys():
						taxonomy=dictTaxidToTaxonomy[taxid]
					else:
						taxonomy=get_desired_ranks(taxid, ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'])
						dictTaxidToTaxonomy[taxid]=taxonomy
					#Ignore uncultured and empty taxonomies
					if len([string for string in taxonomy if 'UNCULTURED' in string.upper()])==0 and len([string for string in taxonomy if 'NONE' == string.upper()])==0:
						#Add to readset that was taxonomised
						setReadsAlignAfterFilteringTaxonomised.add(qname)
						addQnameToTaxDict(qname,taxonomy,dictQnameToTaxonomy)
			ct+=1
			if ct%mark==0:
				pbar.update(ct)


diamondOut.close()
del(dictTaxidToTaxonomy)

print('Exporting fasta format reads which did not pass filtering, taxonomising or could not be aligned',flush=True)
#All fasta sequence IDs which are not in the setReadsAlignAfterFilteringTaxonomised set are written to a fasta file
n=nreadsDedup
mark=round(n/10)
ct=0
with open(filenameOut1, 'w+') as out:
	with gzip.open(DedupFastqFile, 'rt') as f:
		with progressbar.ProgressBar(max_value=n) as pbar:
			pbar.update(ct)
			for seq in SeqIO.parse(f, 'fastq'):
				if seq.id not in setReadsAlignAfterFilteringTaxonomised:
					out.write('>'+str(seq.id)+'\n')
					out.write(str(seq.seq)+'\n')
			ct+=1
			if ct%mark==0:
				pbar.update(ct)
f.close()
out.close()

print('Exporting read numbers',flush=True)
#Calculate length of read sets
ndictBS=len(dictBS)
nsetReadsAlignAfterFiltering=len(setReadsAlignAfterFiltering)
nsetReadsAlignAfterFilteringTaxonomised=len(setReadsAlignAfterFilteringTaxonomised)
del(dictBS)
del(setReadsAlignAfterFiltering)
del(setReadsAlignAfterFilteringTaxonomised)
with open(filenameOut2, 'w+') as out:
	out.write('reads_input\treads_trimmed\treads_deduplicated\treads_mapped\treadsMapQualityFiltered\treadsForTaxonomyCalled\tfile\n')
	out.write(str(nreadsAll)+'\t'+str(nreadsTrim)+'\t'+str(nreadsDedup)+'\t'+str(ndictBS)+'\t'+str(nsetReadsAlignAfterFiltering)+'\t'+str(nsetReadsAlignAfterFilteringTaxonomised)+'\t'+str(filePrefix)+'\n')
out.close()

print('Counting taxons',flush=True)
#For all the found taxa in all available ranks count how often they are present (i.e. have reads assigned to them), if they are hit uniquely by a read
dictCts={}
ct=0
n=len(dictQnameToTaxonomy)
mark=round(n/10)
with progressbar.ProgressBar(max_value=n) as pbar:
	for qname, taxonomy in dictQnameToTaxonomy.items():
		for rank, name in taxonomy.items():
			#If a rank contains multiple entries, e.g. two species do assign a readcount to this rank, else continue
			if len(name)==1:
				#Extract the name element as string from the set
				name=str(list(name)[0])
				#As the len(name)==1 also continues for NA values: If a rank is assigned as NA define an empty set - NA counts in a taxonomy therefore will amount to 0
				if name != 'NA':
					if name not in dictCts.keys():
						dictCts[name]=set()
						dictCts[name].add(qname)
					else:
						dictCts[name].add(qname)
				else:
					dictCts[name]=set()
		ct+=1
		if ct%mark==0:
			pbar.update(ct)

#Calculate readcounts by using the length of the dictionary values (i.e. how many qnames) 
dictCts.update({rank:len(reads) for rank, reads in dictCts.items()})

print('Preparing found taxonomies',flush=True)
#Find available taxonomies which are found by at least one read
normSK=findTaxonomies(dictQnameToTaxonomy,['superkingdom'])
normPH=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum'])
normCL=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum','class'])
normOR=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum','class','order'])
normFA=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum','class','order','family'])
normGE=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum','class','order','family','genus'])
normSP=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum','class','order','family','genus','species'])
normST=findTaxonomies(dictQnameToTaxonomy,['superkingdom','phylum','class','order','family','genus','species','strain'])

#Define readnumber columns for output files
suffix=[str(nreadsAll), str(nreadsTrim), str(nreadsDedup), str(ndictBS), str(nsetReadsAlignAfterFiltering), str(nsetReadsAlignAfterFilteringTaxonomised), str(filePrefix)]

print('Exporting count tables: superkingdom',flush=True)
with open(filenameOut3, 'w+') as out1:
	with open(filenameOut11, 'w+') as out2:
		#Define dictionary
		dictionary=normSK
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			#Extract numbers and build output line
			tax=[superkingdom]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normSK)

print('Exporting count tables: phylum',flush=True)
with open(filenameOut4, 'w+') as out1:
	with open(filenameOut12, 'w+') as out2:
		#Define dictionary
		dictionary=normPH
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normPH)

print('Exporting count tables: class',flush=True)
with open(filenameOut5, 'w+') as out1:
	with open(filenameOut13, 'w+') as out2:
		#Define dictionary
		dictionary=normCL
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			clss=list(dictionary[qname]['class'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum, clss]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normCL)

print('Exporting count tables: order',flush=True)
with open(filenameOut6, 'w+') as out1:
	with open(filenameOut14, 'w+') as out2:
		#Define dictionary
		dictionary=normOR
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			clss=list(dictionary[qname]['class'])[0]
			order=list(dictionary[qname]['order'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum, clss, order]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normOR)

print('Exporting count tables: family',flush=True)
with open(filenameOut7, 'w+') as out1:
	with open(filenameOut15, 'w+') as out2:
		#Define dictionary
		dictionary=normFA
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			clss=list(dictionary[qname]['class'])[0]
			order=list(dictionary[qname]['order'])[0]
			family=list(dictionary[qname]['family'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum, clss, order, family]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normFA)

print('Exporting count tables: genus',flush=True)
with open(filenameOut8, 'w+') as out1:
	with open(filenameOut16, 'w+') as out2:
		#Define dictionary
		dictionary=normGE
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			clss=list(dictionary[qname]['class'])[0]
			order=list(dictionary[qname]['order'])[0]
			family=list(dictionary[qname]['family'])[0]
			genus=list(dictionary[qname]['genus'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum, clss, order, family, genus]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normGE)

print('Exporting count tables: species',flush=True)
with open(filenameOut9, 'w+') as out1:
	with open(filenameOut17, 'w+') as out2:
		#Define dictionary
		dictionary=normSP
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			clss=list(dictionary[qname]['class'])[0]
			order=list(dictionary[qname]['order'])[0]
			family=list(dictionary[qname]['family'])[0]
			genus=list(dictionary[qname]['genus'])[0]
			species=list(dictionary[qname]['species'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum, clss, order, family, genus, species]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normSP)

print('Exporting count tables: strain',flush=True)
with open(filenameOut10, 'w+') as out1:
	with open(filenameOut18, 'w+') as out2:
		#Define dictionary
		dictionary=normST
		for qname, taxonomy in dictionary.items():
			#Extract taxons
			superkingdom=list(dictionary[qname]['superkingdom'])[0]
			phylum=list(dictionary[qname]['phylum'])[0]
			clss=list(dictionary[qname]['class'])[0]
			order=list(dictionary[qname]['order'])[0]
			family=list(dictionary[qname]['family'])[0]
			genus=list(dictionary[qname]['genus'])[0]
			species=list(dictionary[qname]['species'])[0]
			strain=list(dictionary[qname]['strain'])[0]
			#Extract numbers and build output line
			tax=[superkingdom, phylum, clss, order, family, genus, species, strain]
			line=tax + [str(dictCts[rank]) for rank in tax] + suffix
			#Output unfiltered
			out1.write('\t'.join(line)+'\n')
			#Output filtered for read threshold
			if float(100*dictCts[tax[-1]]/nsetReadsAlignAfterFilteringTaxonomised) >= LCAreadSupportThreshold:
				out2.write('\t'.join(line)+'\n')
out1.close()
out2.close()
del(normST)

with open("successful_%s.txt"%filePrefix, 'w+') as out:
	out.write(str(filePrefix)+'\n')
	out.write("%s seconds elapsed"%(time.time()-start_time))
out.close()

print('Finished',flush=True)
print("--- %s seconds elapsed ---"%(time.time()-start_time),flush=True)