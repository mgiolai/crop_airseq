library(tidyverse)
library(RCurl)
library(taxonomizr)

#######################################
##DEFINE LOCATIONS
#######################################
setwd("C:/Users/Michael/Desktop/nhm_analysis/airseq/phibase_dl")
taxDB <- "C:/Users/Michael/Desktop/nhm_analysis/airseq/accessionTaxa.sql"

#######################################
##READ IN PHIBASE FILE TO OBTAIN THE REFERENCE INFORMATION FOR GENOME DOWNLOAD
#######################################

#Download PHIbase file manually first
dfPhibase <- read.delim("phi-base_v4-12_2021-09-02.csv",sep=",", fill = T)
dfPhibase$Pathogen.ID <- as.numeric(dfPhibase$Pathogen.ID)

#Taxonomise with taxonomizr phylogeny
dfTax <- data.frame(getTaxonomy(unique(dfPhibase$Pathogen.ID),taxDB))
dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
dfPhibase <- merge(dfPhibase,dfTax,by.x="Pathogen.ID",by.y="taxid",all=T)
write.csv(dfPhibase,"phibase.taxonomised.csv")


#######################################
##DOWNLOAD REFSEQ AND GENBANK TABLES
#######################################

#Download refseq and genbank tables
#download.file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt","assembly_summary_refseq.txt")
#download.file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt","assembly_summary_genbank.txt")

#Read in refseq dataframe
refseq <- read.delim("assembly_summary_refseq.txt",header=F,sep="\t",skip=2,quote="")
colnames(refseq) <- c("assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date")

#Read in genbank dataframe
genbank <- read.delim("assembly_summary_genbank.txt",header=F,sep="\t",skip=2,quote="")
colnames(genbank) <- c("assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date")

#Remove entries from genbank that are present in refseq to avoid redundancies
genbank <- genbank %>% dplyr::filter(!assembly_accession %in% refseq$gbrs_paired_asm)

#Merge the assembly databases
df <- rbind(refseq,genbank)

#Taxonomise with taxonomizr phylogeny
dfTax <- data.frame(getTaxonomy(unique(df$species_taxid),taxDB))
dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
df <- merge(df,dfTax,by.x="species_taxid",by.y="taxid",all=T)

#######################################
##FILTER THE TABLES
#######################################

#Filter for species present among the PHIbase pathogens and add P. triticina
df <- df %>%
  dplyr::filter(species %in% c(dfPhibase$species,"Puccinia triticina"))

#Only take rows with an https entry
df <- df[grep("https://",df$ftp_path),]

#Export the filtered and taxonomised refseq genbank dataframe
write.table(df,"refseq.genbank.taxonomised.txt",sep="\t",row.names=F,quote=F)

#Query if there is exactly one reference genome, if not, determine the longest
dfAccessions <- c()
for(sp in setdiff(sort(unique(df$species)),dfAccessions$species)){
  print(sp)
  #Subset for species and reference genomes
  x <- df %>% dplyr::filter(species==sp) %>% dplyr::filter(refseq_category=="reference genome")
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Reference Genome"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Reference Genome"))
  }
}

#Query if there is exactly one representative genome, if not, determine the longest
for(sp in setdiff(sort(unique(df$species)),dfAccessions$species)){
  print(sp)
  #Subset for species and representative genomes
  x <- df %>% dplyr::filter(species==sp) %>% dplyr::filter(refseq_category=="representative genome")
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Representative Genome"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Representative Genome"))
  }
}

#Query if there is exactly one complete genome, if not, determine the longest
for(sp in setdiff(sort(unique(df$species)),dfAccessions$species)){
  print(sp)
  #Subset for species and complete genomes
  x <- df %>% dplyr::filter(species==sp) %>% dplyr::filter(assembly_level=="Complete Genome")
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Complete Genome"))
  #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Complete Genome"))
  }
}

#Query if there is exactly one genome, if not, determine the longest
for(sp in setdiff(sort(unique(df$species)),dfAccessions$species)){
  print(sp)
  #Subset for species
  x <- df %>% dplyr::filter(species==sp)
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Other"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    dfAccessions <- rbind(dfAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Other"))
  }
}
write.csv(dfAccessions,"downloaded_accessions.csv")

#######################################
##CREATE A DOWNLOAD FILE
#######################################

df <- df %>% dplyr::filter(assembly_accession %in% dfAccessions$accession)
write.csv(df,"downloaded_refseq_genbank.csv")

#Add download links
df$dlLink <- paste0(df$ftp_path,"/",paste0(sapply(strsplit(df$ftp_path,"/",fixed=T),tail,1),"_genomic.fna.gz"))
df$dlLinkMD5 <- gsub("ftp:","https:",paste0(df$ftp_path,"/md5checksums.txt"))

#Export download file
write.table(df %>% dplyr::select(species,species_taxid,assembly_accession,dlLinkMD5,dlLink),file="download.txt",sep="\t",row.names=F,col.names=F,quote=F)
