library(tidyverse)
library(ggsci)
library(ggvenn)
library(reshape2)
library(treemap)
library(taxonomizr)
library(viridis)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(seqinr)
library(ape)
library(apcluster)
library(pheatmap)
library(data.table)

#######################################
##PREPARE TAXONOMIZR DATABASE
#######################################
setwd("C:/Users/Michael/Desktop")
prepareDatabase("C:/Users/Michael/Desktop/ncbi_taxonomy_db/accessionTaxa.sql")

#######################################
##DEFINE DIRECTORIES
#######################################
parDir <- "C:/Users/Michael/Desktop/nhm_analysis"

# Define other directories
## Directory holding all data
dataDir <- file.path(parDir,"data")
## Directory holding count files
countsDir <- file.path(dataDir,"all-reads_diamond_vs_nr")
## Directory with analysis output folders
analysisDir <- file.path(parDir,"analysis")

########################################
#SPECIFY SAMPLES TO BE ANALYSED
########################################

#Read in sampling datasheet for field experiment
dfSampling <- read.csv(file.path(dataDir,"air-seq_sampling-data_for-R.csv"))

#Read in lab control for background subtraction
labCtrl <- dfSampling %>% dplyr::filter(experiment=="lab.control") %>% pull(sample)
labCtrl <- read.table(file.path(countsDir,paste0("002_6_",labCtrl,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")),sep="\t",header=F)
colnames(labCtrl) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample")

#Read in field samples with increasing sampling times for correlation analysis
dfSamplingSaturation <- dfSampling %>% dplyr::filter(experiment=="field.sampling.altering.samplingtime")
samplesSaturation <- dfSamplingSaturation %>% pull(sample)

#Read in field samples where 60 minute sampling was performed for metagenome analysis
dfSamplingAnalysis <- dfSampling %>% dplyr::filter(collection.minutes==60)
samplesAnalysis <- dfSampling %>% dplyr::filter(collection.minutes==60) %>% pull(sample)

#Read in sampling datasheet for windtunnel sampling
dfWindtunnelSampling <- read.csv(file.path(dataDir,"air-seq_windtunnel-data_for-R.csv"))

#Calculate windunntel harvested spore concentartion
dfWindtunnelSampling$spores <- dfWindtunnelSampling$spore.concentration..spores.ml.* dfWindtunnelSampling$collection.period..minutes. * dfWindtunnelSampling$rate.of.release..ml.minute.
  
#Select 5 metre collection distance samples 
dfWindtunnelSamplingConcSeries <- dfWindtunnelSampling %>% dplyr::filter(distance.from.source..m.=="5" & collection.period..minutes.=="10") %>% dplyr::filter(!rate.of.release..ml.minute.=="0") %>% dplyr::mutate(mio.spores=spores/1e6)
dfWindtunnelSampling10m <- dfWindtunnelSampling %>% dplyr::filter(distance.from.source..m.=="10")

##Load weather data
dfWeather <- read.csv(file.path(dataDir,"weather_data.csv"), header=T)

##Reformat time column
dfWeather$Time <- as.POSIXct(dfWeather$Time, format = "%d/%m/%Y %H:%M", tz="Europe/London")

##Add date column
dfWeather$date <- as.Date(dfWeather$Time)


########################################
#ANALYSE LOAD OF GENUS HOMO IN SAMPLES
########################################

#Load counts
dfTaxCalls <- data.frame()
for(s in dfSampling$sample){
  file <- paste0("002_6_",s,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    print(file)
    #Extract collection minutes time
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
for(s in dfWindtunnelSamplingConcSeries$sample){
  file <- paste0("002_6_airseq_wt_sample-",s,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    print(file)
    #Extract collection minutes time
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample")
dfTaxCalls <- dfTaxCalls %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(nreads.ge.norm=nreads.ge/(nreads.deduplicated/1e6)) %>%
  dplyr::ungroup()

#Visualise unfiltered data
treemap(dfTaxCalls %>%
          dplyr::select(Order,Genus,nreads.ge.norm,sample) %>%
          unique() %>%
          dplyr::group_by(Genus) %>%
          dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
          dplyr::select(Order,Genus,sumCts) %>%
          unique(),
        index = c("Order","Genus"), vSize = "sumCts", title = "",fontsize.labels = 20,fontface.labels=4)

#Calculate human load unfiltered data per sample
tmp <- dfTaxCalls %>%
  dplyr::select(Order,Genus,nreads.ge.norm,sample) %>%
  unique() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
  dplyr::select(sample, Genus, sumCts) %>%
  dplyr::filter(Genus == "Homo") %>%
  unique()

########################################
#LOAD AND ANALYSE WINDTUNNEL SAMPLES
########################################
dir.create(file.path(analysisDir,"01_windtunnel_genera"))
setwd(file.path(analysisDir,"01_windtunnel_genera"))

#Load counts
dfTaxCalls <- data.frame()
for(s in dfWindtunnelSamplingConcSeries$sample){
  file <- paste0("002_6_airseq_wt_sample-",s,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    print(file)
    t <- dfWindtunnelSampling %>% dplyr::filter(sample==s) %>% pull(spores)
    #Extract collection minutes time
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    df$spores <- t
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample","spores")

#Normalise reads RPM
dfTaxCalls <- dfTaxCalls %>% dplyr::mutate(nreads.ge.norm=nreads.ge/(nreads.deduplicated/1e6))
write.csv(dfTaxCalls,"dfTaxCalls.csv",row.names=F)

#Table with readnumbers and released spores
write.csv(dfTaxCalls %>%
            dplyr::select(sample,nreads.sequenced,nreads.trimmed,nreads.deduplicated,nreads.mapped,nreads.filtered,nreads.taxonomised,spores)%>%
            unique(),"readnumbers.csv")

#Plot the number of reads for each collection time-point from sequenced to taxonomised
pdf("readnumbers_sporenumbers.pdf")
tmp <- dfTaxCalls %>%
  dplyr::select(sample,nreads.sequenced,nreads.trimmed,nreads.deduplicated,nreads.mapped,nreads.filtered,nreads.taxonomised,spores) %>%
  unique() %>%
  group_by(sample) %>%
  dplyr::mutate(Sequenced=sum(nreads.sequenced)) %>%
  dplyr::mutate(Trimmed=sum(nreads.trimmed)) %>%
  dplyr::mutate(Deduplicated=sum(nreads.deduplicated)) %>%
  dplyr::mutate(Aligned=sum(nreads.mapped)) %>%
  dplyr::mutate(Filtered=sum(nreads.filtered)) %>%
  dplyr::mutate(Taxonomised=sum(nreads.taxonomised)) %>%
  dplyr::select(sample,Sequenced,Trimmed,Deduplicated,Aligned,Filtered,Taxonomised,spores) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(spores=spores/1e6) %>%
  reshape2::melt(.,id.var=c("sample","spores")) %>%
  unique()

ggplot() + 
  geom_boxplot(data=tmp,aes(x=variable, y=value/1e6,fill=variable),width=0.5) +
  geom_point(data=tmp,aes(x=variable, y=value/1e6),shape=21,size=3) +
  geom_line(data=tmp,aes(x=variable, y=value/1e6,group=sample),linetype="dashed",alpha=0.5) +
  geom_text_repel(data=tmp %>% dplyr::filter(variable=="Sequenced"),aes(x=variable, y=value/1e6,group=sample,label=spores, fontface=2),
                  force             = 2,
                  nudge_x           = -0.50,
                  direction         = "y",
                  hjust             = 0,
                  segment.size      = 0.2,
                  segment.curvature = -0.1) +
  scale_fill_npg() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(fill = "") +
  xlab("") +
  ylab ("Reads [M]") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()


#Treemap of all genera detected in the windtunnel
pdf("windtunnel_treemap_order_genus_background-only.pdf")
treemap(dfTaxCalls %>%
          dplyr::filter(spores==0) %>%
          dplyr::select(Order,Genus,nreads.ge.norm,sample) %>%
          unique() %>%
          dplyr::group_by(Genus) %>%
          dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
          dplyr::select(Order,Genus,sumCts) %>%
          unique(),
        index = c("Order","Genus"), vSize = "sumCts", title = "",fontsize.labels = 20,fontface.labels=4)
treemap(dfTaxCalls %>%
          dplyr::filter(spores==0) %>%
          dplyr::select(Genus,nreads.ge.norm,sample) %>%
          unique() %>%
          dplyr::group_by(Genus) %>%
          dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
          dplyr::select(Genus,sumCts) %>%
          unique(),
        index = c("Genus"), vSize = "sumCts", title = "",fontsize.labels = 20,fontface.labels=4)
dev.off()

write.csv(dfTaxCalls %>%
            dplyr::filter(spores==0) %>%
            dplyr::select(Order,Genus,nreads.ge.norm,sample) %>%
            unique() %>%
            dplyr::group_by(Genus) %>%
            dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
            dplyr::select(Genus,sumCts) %>%
            unique(),"windtunnel_treemap_order_genus_background-only.csv")

pdf("windtunnel_treemap_order_genus_with-background.pdf")
treemap(dfTaxCalls %>%
          dplyr::select(Order,Genus,nreads.ge.norm,sample) %>%
          unique() %>%
          dplyr::group_by(Genus) %>%
          dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
          dplyr::select(Order,Genus,sumCts) %>%
          unique(),
        index = c("Order","Genus"), vSize = "sumCts", title = "",fontsize.labels = 20,fontface.labels=4)
treemap(dfTaxCalls %>%
          dplyr::select(Genus,nreads.ge.norm,sample) %>%
          unique() %>%
          dplyr::group_by(Genus) %>%
          dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
          dplyr::select(Genus,sumCts) %>%
          unique(),
        index = c("Genus"), vSize = "sumCts", title = "",fontsize.labels = 20,fontface.labels=4)
dev.off()

#Load the genera in the 0 spore sample (i.e. the windtunnel background)
airTunnelBackground <- dfTaxCalls %>% dplyr::filter(spores==0) %>% pull(Genus)

#Show treemaps for each spore concentration with removed background
for(sp in sort(setdiff(unique(dfTaxCalls$spores),0))){
  pdf(paste0("windtunnel_treemap_order_genus_",sp,"spores_without-background.pdf"))
  treemap(dfTaxCalls %>%
            dplyr::filter(spores==sp) %>%
            dplyr::filter(!Genus %in% airTunnelBackground) %>%
            dplyr::select(Order,Genus,nreads.ge.norm,sample) %>%
            unique() %>%
            dplyr::group_by(Genus) %>%
            dplyr::mutate(sumCts=sum(nreads.ge.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(sumCts=100*sumCts/sum(nreads.ge.norm)) %>%
            dplyr::select(Order,Genus,sumCts) %>%
            unique(),
          index = c("Order","Genus"),vSize="sumCts",title=sp,fontsize.labels = 30,fontface.labels=4)
  dev.off()
}

########################################
#LOAD WINDTUNNEL SAMPLES AND PLOT AT WHICH CONCENTRATION BACILLUS ARE FOUND
########################################
dir.create(file.path(analysisDir,"02_windtunnel_genera-quantitation"))
setwd(file.path(analysisDir,"02_windtunnel_genera-quantitation"))

#Create a barplot showing the saturation profile of samples detecting Bacillus reads,
#to split up the samples into separate bars per concentration, assign ascending numbers to samples per spore concentration
tmp <- data.frame()
for(sp in sort(unique(dfTaxCalls$spores))){
  ct=1
  samples <- dfTaxCalls %>% dplyr::filter(spores==sp) %>% pull(sample) %>% unique()
  for(sm in samples){
    tmp <- rbind(tmp,dfTaxCalls %>% dplyr::filter(sample==sm) %>% dplyr::mutate(n=ct) %>% dplyr::mutate(spores=spores/1e6))
    ct <- ct +1
  }
}

#Correlation analysis between spores and percent-reads of Bacillus in samples
samples <- tmp %>% 
  dplyr::select(sample,spores,n) %>%
  dplyr::mutate(percentreads=0) %>%
  unique()
corAnalysis <- tmp %>%
  dplyr::select(sample,spores,nreads.ge.norm,Genus,n) %>%
  unique() %>%
  #Sum of samples with equal spore release amounts
  dplyr::group_by(sample,spores) %>%
  dplyr::mutate(nreads.ge.norm.sum=sum(nreads.ge.norm)) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample,spores,nreads.ge.norm,nreads.ge.norm.sum,Genus,n) %>%
  #Percentage reads per spore release concentration
  dplyr::group_by(sample,Genus,spores) %>%
  dplyr::mutate(percentreads=100*sum(nreads.ge.norm)/nreads.ge.norm.sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::filter(Genus=="Bacillus") %>%
  dplyr::select(sample,spores,n,percentreads) %>%
  unique()
corAnalysis <- rbind(corAnalysis,samples %>% dplyr::filter(!sample %in% corAnalysis$sample))
round(cor(corAnalysis$spores, corAnalysis$percentreads,method = c("pearson")),2)
round(cor(corAnalysis$spores, corAnalysis$percentreads,method = c("pearson"))^2,2)
rm(corAnalysis)

#Plot RPM vs % Bacillus 
samples <- tmp %>% 
  dplyr::select(sample,spores,n) %>%
  dplyr::mutate(percentreads=0) %>%
  unique()
corAnalysis <- tmp %>%
  dplyr::select(sample,spores,nreads.ge.norm,Genus,n) %>%
  unique() %>%
  #Sum of samples with equal spore release amounts
  dplyr::group_by(sample,spores) %>%
  dplyr::mutate(nreads.ge.norm.sum=sum(nreads.ge.norm)) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample,spores,nreads.ge.norm,nreads.ge.norm.sum,Genus,n) %>%
  #Percentage reads per spore release concentration
  dplyr::group_by(sample,Genus,spores) %>%
  dplyr::mutate(percentreads=100*sum(nreads.ge.norm)/nreads.ge.norm.sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::filter(Genus=="Bacillus") %>%
  dplyr::select(sample,spores,n,percentreads,nreads.ge.norm) %>%
  unique()
plot(corAnalysis$percentreads~corAnalysis$nreads.ge.norm)


pdf("windtunnel_concentration-bars_sample-saturation-by-spores.pdf",h=5,w=6)
samples <- tmp %>% 
  dplyr::select(sample,spores,n) %>%
  dplyr::mutate(percentreads=0) %>%
  unique()
tmp <- tmp %>%
  dplyr::select(sample,spores,nreads.ge.norm,Genus,n) %>%
  unique() %>%
  #Sum of samples with equal spore release amounts
  dplyr::group_by(sample,spores) %>%
  dplyr::mutate(nreads.ge.norm.sum=sum(nreads.ge.norm)) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample,spores,nreads.ge.norm,nreads.ge.norm.sum,Genus,n) %>%
  #Percentage reads per spore release concentration
  dplyr::group_by(sample,Genus,spores) %>%
  dplyr::mutate(percentreads=100*sum(nreads.ge.norm)/nreads.ge.norm.sum) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::filter(Genus=="Bacillus") %>%
  dplyr::select(sample,spores,n,percentreads) %>%
  unique()
tmp <- rbind(tmp,samples %>% dplyr::filter(!sample %in% tmp$sample))

ggplot() + 
  geom_bar(data=tmp,aes(x=as.character(spores), y=percentreads, group=as.character(n), fill=as.character(spores)), position=position_dodge(), stat="identity", color="black") +
  geom_point(data=tmp %>% dplyr::group_by(spores) %>% dplyr::mutate(mean=mean(percentreads)) %>% dplyr::ungroup(),aes(x=as.character(spores), y=mean),fill="darkred",color="black",shape=22,size=3) +
  geom_line(data=tmp %>% dplyr::group_by(spores) %>% dplyr::mutate(mean=mean(percentreads)) %>% dplyr::ungroup(),aes(x=as.character(spores), y=mean, group=n),linetype="dashed",alpha=0.5) +
  theme_bw(base_size = 14) +
  scale_fill_npg() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16)) +
  labs(fill = "Spores [M]") +
  #facet_wrap(~Genus,scales="free_x",nrow=1) +
  xlab("Spores [M]") +
  ylab(expression("[%] "~italic(Bacillus))) +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Percentage table Bacillus loads concentration series
write.csv(tmp %>% dplyr::group_by(spores) %>% dplyr::mutate(mean=mean(percentreads)) %>% dplyr::ungroup(),"bacillus_loads_concentration_series.csv")

########################################
#ARE BACILLI ALSO FOUND AT 10METRE DISTANCE
########################################
#Load counts
dfTaxCalls <- data.frame()
for(s in dfWindtunnelSampling10m$sample){
  file <- paste0("002_6_airseq_wt_sample-",s,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    print(file)
    t <- dfWindtunnelSampling %>% dplyr::filter(sample==s) %>% pull(spores)
    #Extract collection minutes time
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    df$spores <- t
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample","spores")

#Normalise reads RPM
dfTaxCalls <- dfTaxCalls %>% dplyr::mutate(nreads.ge.norm=nreads.ge/(nreads.deduplicated/1e6))  %>% dplyr::filter(!Genus %in% c(labCtrl$Genus,"Homo"))

write.csv(dfTaxCalls %>%
            dplyr::mutate(spores=spores/1e6) %>%
            dplyr::select(sample,spores,nreads.ge.norm,Genus) %>%
            unique() %>%
            #Sum of samples with equal spore release amounts
            dplyr::group_by(sample,spores) %>%
            dplyr::mutate(nreads.ge.norm.sum=sum(nreads.ge.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(sample,spores,nreads.ge.norm,nreads.ge.norm.sum,Genus) %>%
            #Percentage reads per spore release concentration
            dplyr::group_by(sample,Genus,spores) %>%
            dplyr::mutate(percentreads=100*sum(nreads.ge.norm)/nreads.ge.norm.sum) %>%
            dplyr::ungroup() %>%
            unique() %>%
            dplyr::filter(Genus=="Bacillus") %>%
            dplyr::select(sample,spores,percentreads) %>%
            unique(),"windtunnel_10m-sample_percentages.csv")


########################################
#LOAD FIELD TIME SERIES SAMPLING DATA
########################################
dir.create(file.path(analysisDir,"03_field-sampling-times"))
setwd(file.path(analysisDir,"03_field-sampling-times"))

#Load counts
dfTaxCalls <- data.frame()
for(s in samplesSaturation){
  file <- paste0("002_6_",s,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    #Extract collection minutes time
    t <- dfSampling %>% dplyr::filter(sample==s) %>% pull(collection.minutes)
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    df$collection.minutes <- t
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample","collection.minutes")

#Remove genera which have been found in the lab control
dfTaxCalls <- dfTaxCalls %>% dplyr::filter(!Genus %in% c(labCtrl$Genus,"Homo"))

#Table with readnumbers
write.csv(dfTaxCalls %>%
            dplyr::select(sample,nreads.sequenced,nreads.trimmed,nreads.deduplicated,nreads.mapped,nreads.filtered,nreads.taxonomised,collection.minutes) %>% unique(),"readnumbers.csv")

#Normalise to RPM
dfTaxCalls <- dfTaxCalls %>%
  dplyr::group_by(collection.minutes,Genus) %>%
  dplyr::mutate(nreads.ge.norm=sum(nreads.ge)/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup()

#Number of genera per time-point
write.csv(dfTaxCalls %>% dplyr::select(Genus,sample,collection.minutes) %>% unique() %>% dplyr::group_by(sample) %>% dplyr::mutate(Genus=length(unique(Genus))) %>% dplyr::ungroup() %>% unique(),
          "found-genus.numbers.csv")

#Raw TaxCalls export
write.csv(dfTaxCalls,"dfTaxCalls.csv",row.names=F)

pdf("time-increase-samples_readnumbers.pdf")
tmp <- dfTaxCalls %>%
  dplyr::select(sample,nreads.sequenced,nreads.trimmed,nreads.deduplicated,nreads.mapped,nreads.filtered,nreads.taxonomised,collection.minutes) %>%
  unique() %>%
  group_by(sample) %>%
  dplyr::mutate(Sequenced=sum(nreads.sequenced)) %>%
  dplyr::mutate(Trimmed=sum(nreads.trimmed)) %>%
  dplyr::mutate(Deduplicated=sum(nreads.deduplicated)) %>%
  dplyr::mutate(Aligned=sum(nreads.mapped)) %>%
  dplyr::mutate(Filtered=sum(nreads.filtered)) %>%
  dplyr::mutate(Taxonomised=sum(nreads.taxonomised)) %>%
  dplyr::select(sample,Sequenced,Trimmed,Deduplicated,Aligned,Filtered,Taxonomised,collection.minutes) %>%
  dplyr::ungroup() %>%
  reshape2::melt(.,id.var=c("sample","collection.minutes")) %>%
  unique()

ggplot() + 
  geom_boxplot(data=tmp,aes(x=variable, y=value/1e6,fill=variable),width=0.5) +
  geom_point(data=tmp,aes(x=variable, y=value/1e6),shape=21,size=3) +
  geom_line(data=tmp,aes(x=variable, y=value/1e6,group=sample),linetype="dashed",alpha=0.5) +
  geom_text_repel(data=tmp %>% dplyr::filter(variable=="Sequenced"),aes(x=variable, y=value/1e6,group=sample,label=collection.minutes, fontface=2),
                  force             = 2,
                  nudge_x           = -0.50,
                  direction         = "y",
                  hjust             = 0,
                  segment.size      = 0.2,
                  segment.curvature = -0.1) +
  scale_fill_npg() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(fill = "") +
  xlab("") +
  ylab ("Reads [M]") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Genera found at x sampling time
df <- dfTaxCalls %>% dplyr::select(collection.minutes,nreads.ge.norm,Genus)
df <- reshape2::dcast(data = df,formula = Genus~collection.minutes,value.var ="nreads.ge.norm",fun.aggregate=mean)
row.names(df) <- df$Genus
df <- df[seq(2,ncol(df))]
write.csv(df,"table-for-correlation.csv")
pdf("distance-matrix-sampling-time-similarity-analysis.pdf")
cor(df,use = "pairwise.complete.obs") %>%
  reshape2::melt() %>%
  ggplot(aes(x=factor(Var1,levels=c(5,10,30,60,120)),y=factor(Var2,levels=c(5,10,30,60,120)),fill = value)) +
  geom_tile() +
  geom_text(aes(label=as.character(round(value,2))),color="black",size=6) +
  scale_fill_gradient2(low="#E64B35B2", high="#3C5488B2", guide="colorbar",limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=17),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust = 0.5, size=20),
        axis.text.y=element_text(size=20),
        strip.text = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="right") +
  labs(fill = "") +
  xlab("") +
  ylab ("") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

write.csv(dfTaxCalls %>% dplyr::select(collection.minutes,Genus) %>% dplyr::group_by(collection.minutes) %>%
            dplyr::mutate(nGenera=length(unique(Genus))) %>%
            dplyr::ungroup() %>%
            dplyr::select(collection.minutes,nGenera) %>%
            unique(),"number_genera_collection_time.csv")


#Sampling time vs reads vs number of genera
pdf("observed-Genera_vs_sequencing-depth.pdf",h=6,w=6)
dfTaxCalls %>%
  dplyr::select(sample,nreads.sequenced,collection.minutes,Genus) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(nGenera=length(unique(Genus))) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample,nreads.sequenced,collection.minutes,nGenera) %>%
  unique() %>%
  ggplot() +
  geom_point(aes(x=nreads.sequenced/1000000, y=nGenera, fill=as.character(collection.minutes)),size=6,shape=21) +
  geom_text_repel(aes(x=nreads.sequenced/1000000, y=nGenera, label=paste0(collection.minutes," min."), fontface=2, size=4),
                  force             = 2,
                  nudge_x           = -0.30,
                  #direction         = "y",
                  hjust             = 0,
                  segment.size      = 0.2,
                  segment.curvature = -1) +
  scale_fill_npg() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=17),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust = 0.5, size=20),
        axis.text.y=element_text(size=20),
        strip.text = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(fill = "") +
  xlab("Reads sequenced [M]") +
  ylab ("# Genera") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

write.csv(dfTaxCalls %>%
            dplyr::select(collection.minutes,Genus,sample) %>%
            dplyr::group_by(collection.minutes,sample) %>%
            dplyr::mutate(nGenera=length(unique(Genus))) %>%
            dplyr::ungroup() %>%
            dplyr::select(collection.minutes,sample,nGenera) %>%
            unique(),"number_genera_over_time.csv")

#Associate sampling time with library amount
pdf("observed-Genera_vs_ng-sequencing-library.pdf",h=6,w=6)
dfTaxCalls %>%
  dplyr::select(sample,collection.minutes,Genus,nreads.sequenced,nreads.taxonomised) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(nGenera=length(unique(Genus))) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample,collection.minutes,nGenera,nreads.sequenced,nreads.taxonomised) %>%
  unique() %>%
  dplyr::mutate(concentration=c(0.48525,0.38291,0,0.16074,3.9575,3.5156,3.8102,16.092,17.33,18.41)) %>%
  ggplot() +
  geom_point(aes(y=nGenera, x=concentration, fill=as.character(collection.minutes)),size=6,shape=21) +
  geom_smooth(aes(y=nGenera, x=concentration),linetype="dashed",colour="grey34",size=1,se=F,method='loess') +
  geom_text_repel(aes(y=nGenera, x=concentration, label=paste0(collection.minutes," min."), fontface=2, size=4)) +
  scale_fill_npg() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=17),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=90, vjust = 0.5, size=20),
        axis.text.y=element_text(size=20),
        strip.text = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(fill = "") +
  ylab("No. Genera") +
  xlab ("Nextera library [ng]") +
  ggtitle("") +
  geom_vline(xintercept=3.5,linetype="dashed") +
  scale_x_continuous(breaks=c(0,3.5,5,10,15,20)) +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()


#Show Venn diagram overlap of species
pdf("genera-overlaps-times_all.pdf")
#30, 60, 120 minutes
x <- list(A=dfTaxCalls %>% dplyr::filter(collection.minutes==30) %>% pull(Genus) %>% unique(),
          B=dfTaxCalls %>% dplyr::filter(collection.minutes==60) %>% pull(Genus) %>% unique(),
          C=dfTaxCalls %>% dplyr::filter(collection.minutes==120) %>% pull(Genus) %>% unique())
names(x) <- c("30 min","60 min","120 min")
ggvenn(x,
       fill_color = c("#4DBBD5B2","#EFC000FF","#DC0000B2"),
       stroke_size = 0.5,
       set_name_size = 4) -> p
plot(p)

#10, 30, 60, 120 minutes
x <- list(A=dfTaxCalls %>% dplyr::filter(collection.minutes==10) %>% pull(Genus) %>% unique(),
          B=dfTaxCalls %>% dplyr::filter(collection.minutes==30) %>% pull(Genus) %>% unique(),
          C=dfTaxCalls %>% dplyr::filter(collection.minutes==60) %>% pull(Genus) %>% unique(),
          D=dfTaxCalls %>% dplyr::filter(collection.minutes==120) %>% pull(Genus) %>% unique())
names(x) <- c("10 min","30 min","60 min","120 min")
ggvenn(x,
       fill_color = c("#4DBBD5B2","#EFC000FF","#DC0000B2","#3C5488B2"),
       stroke_size = 0.5,
       set_name_size = 4) -> p
plot(p)
dev.off()

#Plot number of reads in samples
pdf("reads_taxonomised_time-point.pdf")
dfTaxCalls %>%
  dplyr::select(collection.minutes,sample,nreads.taxonomised) %>%
  unique() %>%
  dplyr::group_by(collection.minutes) %>%
  dplyr::mutate(sumReads=sum(nreads.taxonomised)/1e6) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x=factor(collection.minutes,levels=c(5,10,30,60,120)), y=sumReads, fill=as.character(collection.minutes))) + 
  geom_bar(stat="identity", position="identity") +
  theme_bw(base_size = 14) +
  scale_fill_npg() +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(fill = "") +
  xlab("Collection time") +
  ylab ("Reads taxonomised [M]") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

write.csv(dfTaxCalls %>%
            dplyr::select(collection.minutes,sample,nreads.taxonomised) %>%
            unique() %>%
            dplyr::group_by(collection.minutes) %>%
            dplyr::mutate(sumReads=sum(nreads.taxonomised)) %>%
            dplyr::ungroup(),"reads_taxonomised_timepoint.csv")


########################################
#ANALYSE FIELD TIMECOURSE
########################################
dir.create(file.path(analysisDir,"04_field-sampling-timecourse_presence"))
setwd(file.path(analysisDir,"04_field-sampling-timecourse_presence"))

#Load counts
dfTaxCalls <- data.frame()
for(sample in samplesAnalysis){
  file <- paste0("002_6_",sample,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample")

#Associate samples with date and time
dfTaxCalls <- merge(dfTaxCalls,dfSamplingAnalysis[c("sample","date")],by="sample")
dfTaxCalls$date <-  as.Date(dfTaxCalls$date, format = "%d/%m/%Y")

#Normalise to RPM per replicate
tmp <- dfTaxCalls %>%
  dplyr::group_by(sample,Genus) %>%
  dplyr::mutate(nreads.ge=sum(nreads.ge)) %>%
  dplyr::mutate(nreads.ge.norm=nreads.ge/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(!c(nreads.deduplicated)) %>%
  unique() %>%
  dplyr::select(!nreads.ge) %>%
  dplyr::select(sample,date,Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm) %>% 
  unique()
write.csv(tmp,"replicate_date_RPM_per_genu.csv")

#Normalise to RPM per date
dfTaxCalls <- dfTaxCalls %>%
  dplyr::group_by(date,Genus) %>%
  dplyr::mutate(nreads.ge.norm=sum(nreads.ge)/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup()

#Remove genera which have been found in the lab control
dfTaxCalls <- dfTaxCalls %>% dplyr::filter(!Genus %in% c(labCtrl$Genus,"Homo"))

#TaxCalls raw export
write.csv(dfTaxCalls,"dfTaxCalls.csv",row.names=F)

#Table with readnumbers
write.csv(dfTaxCalls %>%
            dplyr::select(sample,nreads.sequenced,nreads.trimmed,nreads.deduplicated,nreads.mapped,nreads.filtered,nreads.taxonomised) %>% unique(),"readnumbers.csv")

pdf("readnumbers-field-time-series.pdf")
#Plot number of reads in samples
dfTaxCalls %>%
  dplyr::select(date,nreads.sequenced,nreads.trimmed,nreads.deduplicated,nreads.mapped,nreads.filtered,nreads.taxonomised) %>%
  unique() %>%
  group_by(date) %>%
  dplyr::mutate(Sequenced=sum(nreads.sequenced)) %>%
  dplyr::mutate(Trimmed=sum(nreads.trimmed)) %>%
  dplyr::mutate(Deduplicated=sum(nreads.deduplicated)) %>%
  dplyr::mutate(Aligned=sum(nreads.mapped)) %>%
  dplyr::mutate(Filtered=sum(nreads.filtered)) %>%
  dplyr::mutate(Taxonomised=sum(nreads.taxonomised)) %>%
  dplyr::select(date,Sequenced,Trimmed,Deduplicated,Aligned,Filtered,Taxonomised) %>%
  dplyr::ungroup() %>%
  reshape2::melt(.,id.var=c("date")) %>%
  unique() %>%
  ggplot() + 
  geom_boxplot(aes(x=variable, y=value/1e6,fill=variable)) +
  geom_point(aes(x=variable, y=value/1e6),shape=21,size=3) +
  geom_line(aes(x=variable, y=value/1e6,group=date),linetype="dashed",alpha=0.5) +
  scale_fill_npg() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(fill = "") +
  xlab("") +
  ylab ("Reads [M]") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Percentages of Genus taxonomic levels found over entire series
write.csv(dfTaxCalls %>%
            dplyr::select(Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm,sample) %>%
            unique() %>%
            dplyr::group_by(Superkingdom,Phylum,Class,Order,Family,Genus,) %>%
            dplyr::mutate(percent=sum(nreads.ge.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm,percent) %>%
            dplyr::mutate(percent=100*percent/sum(nreads.ge.norm)) %>%
            unique() %>%
            dplyr::select(Superkingdom,Phylum,Class,Order,Family,Genus,percent) %>%
            unique(),
          "percent-entire-series-genus.csv")

#Treemaps of most abundant Phyla, Order, Species
pdf("jic-time-series-treemap_phylum_phylum_genus.pdf")
treemap(dfTaxCalls %>%
          dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
          unique() %>%
          dplyr::group_by(Phylum,Genus) %>%
          dplyr::mutate(percentReads=sum(nreads.ge.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(percentReads=100*percentReads/sum(nreads.ge.norm)) %>%
          dplyr::select(Phylum,Genus,percentReads) %>%
          unique(), index = c("Phylum","Genus"), vSize = "percentReads", title = "",fontsize.labels = 25,fontface.labels=4)
dev.off()

write.csv(dfTaxCalls %>%
            dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
            unique() %>%
            dplyr::group_by(Phylum,Genus) %>%
            dplyr::mutate(percentReads=sum(nreads.ge.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(percentReads=100*percentReads/sum(nreads.ge.norm)) %>%
            dplyr::select(Phylum,Genus,percentReads) %>%
            unique(),"treemap_table.csv")

#Number of Genera per Phylum
write.csv(dfTaxCalls %>%
            dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
            unique() %>%
            group_by(Phylum) %>%
            dplyr::mutate(nGenera=length(unique(Genus))) %>%
            dplyr::mutate(sumReads=sum(nreads.ge.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(Phylum,nGenera,sumReads) %>%
            unique(),"number-genera-per-phylum.csv")

########################################
#ANALYSE FIELD TIMECOURSE: CORRELATION BETWEEN REPLICATES ON THE SAME DAY
########################################
dir.create(file.path(analysisDir,"05_field-sampling-timecourse_replicate-comparison"))
setwd(file.path(analysisDir,"05_field-sampling-timecourse_replicate-comparison"))

#Load counts
dfTaxCalls <- data.frame()
for(sample in samplesAnalysis){
  file <- paste0("002_6_",sample,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample")

#Associate samples with date and time
dfTaxCalls <- merge(dfTaxCalls,dfSamplingAnalysis[c("sample","date")],by="sample")
dfTaxCalls$date <-  as.Date(dfTaxCalls$date, format = "%d/%m/%Y")

#Normalise to RPM
dfTaxCalls <- dfTaxCalls %>%
  dplyr::group_by(sample,Genus) %>%
  dplyr::mutate(nreads.ge.norm=sum(nreads.ge)/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup()

#Remove genera which have been found in the lab control
dfTaxCalls <- dfTaxCalls %>% dplyr::filter(!Genus %in% c(labCtrl$Genus,"Homo"))

#Prepare dataframe
df <- dfTaxCalls %>%
  dplyr::select(sample,Genus,nreads.ge.norm) %>%
  reshape2::dcast(.,formula=Genus~sample,value.var="nreads.ge.norm",fun.aggregate=mean)
row.names(df) <- df$sample
df <- df[seq(2,ncol(df))]

#Perform PCA
tmp <- df
tmp[is.na(tmp)] <- 0
pca_result <- prcomp(t(tmp),scale.=T)

#Scree plot
plot(pca_result,type="l",main="Scree Plot")

#Extract proportion of variance for PC1 and PC2
V1 <- round(as.data.frame(summary(pca_result)$importance)[2,1]*100,2)
V2 <- round(as.data.frame(summary(pca_result)$importance)[2,2]*100,2)

#Plot PCA
pcaData <- as.data.frame(pca_result$x[, 1:2])
pcaData$sample <- row.names(pcaData)
pdf("PCA_replicates-per-date.pdf",h=12,w=12)
pcaData %>%
  merge(.,dfTaxCalls %>% dplyr::select(sample,date) %>% unique() %>% dplyr::group_by(date) %>% dplyr::mutate(ct=row_number()),by="sample") %>%
  dplyr::mutate(Var1=paste0(date,".rep",ct)) %>%
  dplyr::mutate(Var1=gsub("2015-","",Var1)) %>%
  dplyr::select(sample=Var1,PC1,PC2) %>%
  ggplot(aes(x=PC1,y=PC2)) +
  geom_point(aes(fill=sample),shape=21,size=2) +
  geom_text_repel(aes(label=sample),box.padding = 0.5, point.padding = 1, segment.color = "grey", nudge_x = 0.1, nudge_y = 0.1, max.overlaps = Inf) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x=element_text(angle=90, vjust = 0.5, size=12),
        axis.text.y=element_text(size=12),
        strip.text = element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position="") +
  labs(fill = "") +
  xlab(paste0("PC1 [",V1,"%]")) +
  ylab (paste0("PC2 [",V2,"%]")) +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Perform correlation analysis
tmp <- cor(df,use = "pairwise.complete.obs") %>%
  reshape2::melt() %>%
  merge(.,dfTaxCalls %>% dplyr::select(sample,date) %>% unique() %>% dplyr::group_by(date) %>% dplyr::mutate(ct=row_number()),by.x="Var1",by.y="sample") %>%
  dplyr::mutate(Var1=paste0(date,".rep",ct)) %>%
  dplyr::select(!c(date,ct)) %>%
  merge(.,dfTaxCalls %>% dplyr::select(sample,date) %>% unique() %>% dplyr::group_by(date) %>% dplyr::mutate(ct=row_number()),by.x="Var2",by.y="sample") %>%
  dplyr::mutate(Var2=paste0(date,".rep",ct)) %>%
  dplyr::select(!c(date,ct)) %>%
  dplyr::mutate(Var1=gsub("2015-","",Var1)) %>%
  dplyr::mutate(Var2=gsub("2015-","",Var2)) 

#Plot correlation
pdf("correlation_replicates-per-date.pdf",h=12,w=12)
tmp %>%
  ggplot(aes(x=Var1,Var2,fill = value)) +
  geom_tile() +
  geom_text(aes(label=as.character(round(value,2))),color="black",size=2) +
  scale_fill_gradient2(low="#E64B35B2", high="#3C5488B2", guide="colorbar",limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.x=element_text(angle=90, vjust = 0.5, size=8),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        legend.position="right") +
  labs(fill = "") +
  xlab("") +
  ylab ("") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Analyse range of correlation values for the same day
tmp <- tmp %>% dplyr::filter(!Var1==Var2)
tmp$x1 <- sapply(strsplit(tmp$Var1,".",fixed=T),'[[',1)
tmp$x2 <- sapply(strsplit(tmp$Var2,".",fixed=T),'[[',1)
tmp1 <- tmp %>%
  dplyr::filter(x1==x2) %>%
  dplyr::group_by(x1) %>%
  dplyr::mutate(mean=mean(value)) %>% 
  dplyr::ungroup() %>%
  dplyr::select(x1,mean) %>%
  unique() %>%
  dplyr::mutate(meanofall=mean(mean)) %>%
  dplyr::mutate(medianofall=median(mean))

hist(tmp1$mean)
write.csv(tmp1, "correlation-comparison_samplesof-same-day.csv")

########################################
#PLOT NUMBER OF READS PER PHYLUM OVER TIME
########################################
dir.create(file.path(analysisDir,"06_field-sampling-timecourse_timeplots-genus"))
setwd(file.path(analysisDir,"06_field-sampling-timecourse_timeplots-genus"))

#Load counts
dfTaxCalls <- data.frame()
for(sample in samplesAnalysis){
  file <- paste0("002_6_",sample,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
  if(file %in% list.files(countsDir)){
    df<-read.table(file.path(countsDir,file),header=F,sep='\t')
    dfTaxCalls <- rbind(dfTaxCalls,df)
  }
}
colnames(dfTaxCalls) <- c("Superkingdom","Phylum","Class","Order","Family","Genus","nreads.sk","nreads.phy","nreads.cl","nreads.or","nreads.fa","nreads.ge","nreads.sequenced","nreads.trimmed","nreads.deduplicated","nreads.mapped","nreads.filtered","nreads.taxonomised","sample")

#Remove genera which have been found in the lab control
dfTaxCalls <- dfTaxCalls %>% dplyr::filter(!Genus %in% c(labCtrl$Genus,"Homo"))

#Associate samples with date and time
dfTaxCalls <- merge(dfTaxCalls,dfSamplingAnalysis[c("sample","date")],by="sample")
dfTaxCalls$date <-  as.Date(dfTaxCalls$date, format = "%d/%m/%Y")

#Normalise to RPM per sample
tmp <- dfTaxCalls %>%
  dplyr::group_by(sample,Genus) %>%
  dplyr::mutate(nreads.ge.norm=sum(nreads.ge)/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm,sample,date) %>%
  unique()
write.csv(tmp,"RPM-per-sample_dates-shown.csv")

#Normalise to RPM
dfTaxCalls <- dfTaxCalls %>%
  dplyr::group_by(date,Genus) %>%
  dplyr::mutate(nreads.ge.norm=sum(nreads.ge)/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup()


#Genera numbers per phylum over time
pdf("jic_per-phylum-genera-over-time.pdf",height = 8,width = 10)
orderPhylum<- dfTaxCalls %>% dplyr::select(Superkingdom,Phylum) %>% unique() %>% dplyr::arrange(Superkingdom) %>% pull(Phylum) %>% rev()
orderDate <- allDatesReformat <- dfTaxCalls %>% dplyr::arrange(date) %>% dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>% pull(dateReformat) %>% unique()
dfTaxCalls %>%
  dplyr::group_by(date,Superkingdom,Phylum) %>% 
  dplyr::mutate(nGenus=length(unique(Genus))) %>%
  dplyr::ungroup() %>%
  dplyr::select(date,Superkingdom,Phylum,nGenus) %>%
  unique() %>%
  dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>%
  ggplot() +
  geom_tile(aes(x=factor(dateReformat,levels=orderDate), y=factor(Phylum,levels=orderPhylum), fill=nGenus)) +
  coord_equal() +
  scale_fill_gradientn(trans="log10", colours =  rev(brewer.pal(11, 'RdYlBu'))) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=20),
        axis.text.y=element_text(face = "italic"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position="bottom") +
  labs(fill = "Genera") +
  xlab("") +
  ylab ("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Reads over time
pdf("jic_per-phylum-sumreads-over-time.pdf",height = 8,width = 10)
orderPhylum<- dfTaxCalls %>% dplyr::select(Superkingdom,Phylum) %>% unique() %>% dplyr::arrange(Superkingdom) %>% pull(Phylum) %>% rev()
orderDate <- allDatesReformat <- dfTaxCalls %>% dplyr::arrange(date) %>% dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>% pull(dateReformat) %>% unique()
dfTaxCalls %>%
  dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
  unique() %>%
  #Sum up normalised reads genus to create phylum counts
  dplyr::group_by(date,Phylum) %>%
  dplyr::mutate(nreads.ge.norm=sum(nreads.ge.norm)) %>%
  dplyr::ungroup() %>%
  dplyr::select(date,Phylum,nreads.ge.norm) %>%
  unique() %>%
  dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>%
  ggplot() +
  geom_tile(aes(x=factor(dateReformat,levels=orderDate), y=factor(Phylum,levels=orderPhylum), fill=nreads.ge.norm)) +
  coord_equal() +
  scale_fill_gradientn(trans="log10", colours =  rev(brewer.pal(11, 'RdYlBu')), breaks=c(100,1000,10000,100000)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=20),
        axis.text.y=element_text(face = "italic"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position="bottom") +
  labs(fill = "RPM dedup.") +
  xlab("") +
  ylab ("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Export a table for All, Arthropoda and ascomycots
tmp <- dfTaxCalls %>%
  dplyr::select(date,Genus,nreads.ge.norm) %>%
  unique() %>%
  reshape2::dcast(.,formula=Genus~date,value.var="nreads.ge.norm",fun.aggregate=mean)
tmp <- dfTaxCalls %>%
  dplyr::select(Superkingdom,Phylum,Class,Order,Family,Genus) %>%
  merge(.,tmp,by="Genus") %>%
  unique()
write.csv(tmp,"RPM_matrix.csv")
write.csv(dfTaxCalls %>% dplyr::select(sample,date,Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm) %>% unique(),
          "arthropoda_genus_by-Sample-by-date.csv")
write.csv(dfTaxCalls %>% dplyr::select(sample,date,Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm) %>% unique() %>% dplyr::filter(Phylum=="Arthropoda"),
          "arthropoda_genus_by-Sample-by-date.csv")
write.csv(dfTaxCalls %>% dplyr::select(sample,date,Superkingdom,Phylum,Class,Order,Family,Genus,nreads.ge.norm) %>% unique() %>% dplyr::filter(Phylum=="Ascomycota"),
          "ascomycota_genus_by-Sample-by-date.csv")

pdf("jic_arthropoda-sumreads-over-time.pdf",height = 15,width = 10)
orderDate <- allDatesReformat <- dfTaxCalls %>% dplyr::arrange(date) %>% dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>% pull(dateReformat) %>% unique()
dfTaxCalls %>%
  dplyr::filter(Phylum=="Arthropoda") %>%
  dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
  unique() %>%
  dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>%
  ggplot() +
  geom_tile(aes(x=factor(dateReformat,levels=orderDate), y=Genus, fill=nreads.ge.norm)) +
  coord_equal() +
  scale_fill_gradientn(trans="log10", colours =  rev(brewer.pal(11, 'RdYlBu')), breaks=c(100,1000,10000,100000)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=20),
        axis.text.y=element_text(face = "italic"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position="bottom") +
  labs(fill = "RPM dedup.") +
  xlab("") +
  ylab ("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()


pdf("jic_ascomycota-sumreads-over-time.pdf",height = 12,width = 10)
orderDate <- allDatesReformat <- dfTaxCalls %>% dplyr::arrange(date) %>% dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>% pull(dateReformat) %>% unique()
dfTaxCalls %>%
  dplyr::filter(Phylum=="Ascomycota") %>%
  dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
  unique() %>%
  dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>%
  ggplot() +
  geom_tile(aes(x=factor(dateReformat,levels=orderDate), y=Genus, fill=nreads.ge.norm)) +
  coord_equal() +
  scale_fill_gradientn(trans="log10", colours =  rev(brewer.pal(11, 'RdYlBu')), breaks=c(100,1000,10000,100000)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=20),
        axis.text.y=element_text(face = "italic"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position="bottom") +
  labs(fill = "RPM dedup.") +
  xlab("") +
  ylab ("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

########################################
#SHOW HOW PHIBASE GENERA VARY OVER TIME
########################################
dir.create(file.path(analysisDir,"07_field-sampling-timecourse_phibase-genera"))
setwd(file.path(analysisDir,"07_field-sampling-timecourse_phibase-genera"))

#Load and taxonomise phibase pathogen and host data
df <- read.csv(file.path(dataDir,"phi-base_v4-12_2021-09-02.csv")) %>% dplyr::select(Pathogen.ID,Host.ID)
colnames(df) <- c("taxid.pathogen","taxid.host")
##Call taxa on phibase data
df$taxid.pathogen <- as.numeric(df$taxid.pathogen) #needs to be numeric
dfTax <- data.frame(getTaxonomy(unique(df$taxid.pathogen),'C:/Users/Michael/Desktop/ncbi_taxonomy_db/accessionTaxa.sql'))
dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
##Control taxid content
setdiff(df$taxid.pathogen,dfTax$taxid)
setdiff(dfTax$taxid,df$taxid.pathogen)
length(intersect(dfTax$taxid,df$taxid.pathogen))
length(unique(df$taxid.pathogen))
##Merge to obtain phylogeny
df <- merge(df,dfTax,by.x="taxid.pathogen",by.y="taxid",all=T)
##Call taxa on phibase data
df$taxid.host <- as.numeric(df$taxid.host) #needs to be numeric
dfTax <- data.frame(getTaxonomy(unique(df$taxid.host),'C:/Users/Michael/Desktop/ncbi_taxonomy_db/accessionTaxa.sql'))
dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
##Control taxid content
setdiff(df$taxid.host,dfTax$taxid)
setdiff(dfTax$taxid,df$taxid.host)
length(intersect(dfTax$taxid,df$taxid.host))
length(unique(df$taxid.host))
##Merge to obtain phylogeny
df <- merge(df,dfTax,by.x="taxid.host",by.y="taxid",all=T)
#Exchange column names .x and .y
colnames(df) <- gsub(".x",".pathogen",colnames(df),fixed=T)
colnames(df) <- gsub(".y",".host",colnames(df),fixed=T)

#Write table which pathogen genera are found per phylum
write.table(dfTaxCalls %>%
  dplyr::select(Phylum,Genus) %>% 
  dplyr::filter(Genus %in% df$genus.pathogen) %>%
  unique() %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(Genus=paste(sort(Genus),collapse=", ")) %>%
  dplyr::ungroup() %>%
  unique(),"which-pathogen_genera_per_phylum.tsv",sep="\t")

#Export list of all genera with pathogen / non-pathogen label
tmp <- dfTaxCalls %>% dplyr::select(Phylum,Genus) %>% unique()
tmp$Pathogen <- sapply(tmp$Genus, function(x) ifelse(x %in% df$genus.pathogen,"p.","n.p."))
write.csv(tmp,"phibase_scoring-genera-per-phlyum.csv")

#Calculate how often a Genus is observed and the average number of RPM, match with the phiabase data to identify pathogens
tmp <- dfTaxCalls %>%
  dplyr::select(date,Phylum,Genus,nreads.ge.norm) %>%
  unique() %>%
  #Average change of entire sample over time
  dplyr::group_by(Genus) %>%
  dplyr::mutate(meanRPM=mean(nreads.ge.norm)) %>%
  dplyr::mutate(days=length(unique(date))) %>%
  dplyr::ungroup() %>%
  dplyr::select(Phylum,Genus,days,meanRPM) %>%
  unique()
tmp$Pathogen <- sapply(tmp$Genus, function(x) ifelse(x %in% df$genus.pathogen,"p.","n.p."))
write.csv(tmp,"meanRPM_daysObserved-pathogens.csv")


#Plots for pathogens
pdf("species_and_phibase_days_vs_meanRPM.pdf",width=9,height=4)
#Filter for Phyla containing phibase hits
tmp %>%
  dplyr::filter(Phylum %in% c(tmp %>% dplyr::filter(Pathogen=="p.") %>% dplyr::pull(Phylum) %>% unique())) %>%
  reshape2::melt(id.vars=c("Phylum","Genus","days","Pathogen")) %>%
  ggplot() +
  geom_point(aes(x=days, y=log10(value), colour=Pathogen, shape=Pathogen),size=3) +
  theme_bw(base_size = 14) +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_aaas() +
  scale_shape_manual(values=c(21, 15))+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=16),
        axis.text.y=element_text(size=16),
        strip.text.x = element_text(size = 9),
        strip.text=element_text(face = "italic"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="right") +
  labs(color  = "Pathogen", shape = "Pathogen") +
  xlab("Days observed") +
  ylab ("log10 mean(RPM  dedup.)") +
  ggtitle("") +
  facet_wrap(~Phylum,nrow=2) -> p
plot(p)
dev.off()

########################################
#SHOW HOW PHIBASE SPECIES VARY OVER TIME: BWA MAPPING CUTOFF ANALYSIS
########################################
dir.create(file.path(analysisDir,"08_field-sampling-timecourse_phibase-species_bwa-cutoff-analysis"))
setwd(file.path(analysisDir,"08_field-sampling-timecourse_phibase-species_bwa-cutoff-analysis"))

#Read in pmatch and pident filtered sam files
samplesFiltering <- c()
for(file in grep("pident",list.files(file.path(dataDir,"bwa_airseq_vs_phibase4-12")),value=T)){
  sample <- strsplit(file,".",fixed=T)[[1]][1]
  if(sample %in% dfSamplingAnalysis$sample){
    samplesFiltering <- unique(c(samplesFiltering,sample))
  }
}
#Convert dfSamplingAnalysis to a data.table
dtSamplingAnalysis <- as.data.table(dfSamplingAnalysis)
#Convert date column in dtSamplingAnalysis to Date format
dtSamplingAnalysis[, date := as.Date(date, format = "%d/%m/%Y")]

#Read genome sizes of referenes for filtering
refSize <- read.table(file.path(dataDir,"reference_genomesizes.txt"),header=F,sep="\t")
colnames(refSize) <- c("taxid","size","assembly")

for(i in 1:length(samplesFiltering)){
  #Sample name
  sample <- samplesFiltering[i]
  
  #Output filename
  rdsName <- paste0(sample,".filtered.rds")
  csvName <- paste0(sample,".filtered.csv")

  if(!file.exists(rdsName) | !file.exists(csvName)){
    #If the output filename does not exist, proceed
    samgz <- grep(paste0(sample,"."),grep("pident",list.files(file.path(dataDir,"bwa_airseq_vs_phibase4-12")),value=T),value=T,fixed=T)
    print(samgz)
    
    #Read sam.gz file
    rL <- read_lines(file(file.path(dataDir,"bwa_airseq_vs_phibase4-12",samgz)))
    
    #Keep all lines but the sam header
    rL <- rL[!startsWith(rL, "@")]
    
    #Remove reads shorter than x bp, calculate number of taxonomised reads, remove all reads mapping with a MAPQ of 0
    queryKeep <- intersect(which(nchar(unlist(lapply(rL, function(x) strsplit(x, "\t")[[1]][10])))>=0), # read length
                           which(lapply(rL, function(x) as.numeric(strsplit(x, "\t")[[1]][5]))>=0)) # mapq
    rL <- rL[c(queryKeep)]
    
    #Count number of taxonomised reads
    all.nreads.taxonomised <- length(unique(unlist(lapply(rL, function(x) strsplit(x, "\t")[[1]][1]))))
    
    #Count number of taxonomised bases
    all.bases.taxonomised <- sum(unlist(lapply(rL, function(x) nchar(strsplit(x, "\t")[[1]][10]))))
    
    #Extract all the taxids that are being hit by a read, if there are multiple for a single read, discard the read
    if(length(rL)>0){
      #Extract all taxids and make them unique to see if length==1 (this will be used for filtering reads that map to more than one taxid)
      taxids <- lapply(str_extract_all(rL, "(?<=taxid\\|)\\d+"), unique)
      
      #Calculate which read hits more than one taxid based on SA:Z and XA:Z flags present in the line as these reads will be removed
      taxidsLen <- lapply(taxids, length)
      #Calculate which taxids are of length 1 as these will be kept
      queryKeep <- which(unlist(taxidsLen)==1)
      #Only keep lines with unique hits to the same taxid
      rL <- rL[c(queryKeep)]
      
      #Only keep taxids of reads where a single taxid was hit - i.e. this removes reads that map to more than one species
      taxids <- as.numeric(unlist(taxids[c(queryKeep)]))
      
      if(length(rL)>0){
        #Create sam data table using query and taxid information
        k <- strsplit(rL, "\t")
        sam <- data.table(query=sapply(k, function(x) as.character(x[1])),
          mapq=sapply(k, function(x) as.numeric(x[5])),
          reference=sapply(k, function(x) as.character(x[3])),
          start=sapply(k, function(x) as.numeric(x[4])),
          orientation=sapply(k, function(x) as.numeric(x[2])),
          cigar=sapply(k, function(x) as.character(x[6])),
          taxid=as.numeric(taxids),
          #readLen=nchar(sapply(k, function(x) as.character(x[10]))),
          all.nreads.taxonomised=all.nreads.taxonomised)
        
        #Calculate which queries will be kept as they only find a single species
        sam <- sam[, .(num_taxid = uniqueN(taxid), query, start, reference, cigar, orientation, all.nreads.taxonomised, all.bases.taxonomised), by = taxid][num_taxid == 1]
        sam[, num_taxid := NULL]
        
        #Calculate number of reads per species
        sam[, nreads.sp := length(unique(query)), by = taxid]
        
        #Calculate percent reads for abundance filtering
        sam[, percentReads := 100*nreads.sp/all.nreads.taxonomised, by = taxid]
        
        #Remove human taxid
        sam <- sam[!taxid==9606]
        
        #Split reverse mapping reads
        tmp <- sam[orientation==16]
        
        #List holding information which references bases are covered
        regionsCovered <- list()
        #Calculate which regions are covered by at least one reverse read
        for(ref in unique(tmp$reference)){
          print(ref)
          positionCovered <- integer(0)
          tmp1 <- tmp[reference==ref,][,cigar,start]
          tmp1[, count := .N, by = .(cigar, start)]
          tmp1 <- unique(tmp1)
          for(row in 1:nrow(tmp1)){
            start <- tmp1[row,]$start+1
            cigarString <- tmp1[row,]$cigar
            count <- tmp1[row,]$count
            #Regular expression to split the CIGAR string into operations
            cigarSplit <- regmatches(cigarString, gregexpr("([0-9]+)([MIDNSHPX=])", cigarString))[[1]]
            numbers <- as.numeric(gsub("[A-Z=]", "", cigarSplit))
            characters <- gsub("[0-9]", "", cigarSplit)
            for(type in 1:length(characters)){
              if(characters[type] %in% c("M","X","=")){
                end=tail(start-1:numbers[type],n=1)
                positionCovered <- c(positionCovered,rep(paste0(start,":",end),count))
                start=end
              }else if(characters[type] %in% c("I")){
                start=start
              }else{
                start=tail(start-1:numbers[type],n=1)              }
            }
          }
          regionsCovered[[ref]] <- positionCovered
        }
          
        #Split forward mapping reads
        tmp <- sam[!orientation==16]
        
        #Calculate which regions are covered by at least one forward read
        for(ref in unique(tmp$reference)){
          print(ref)
          positionCovered <- integer(0)
          tmp1 <- tmp[reference==ref,][,cigar,start]
          tmp1[, count := .N, by = .(cigar, start)]
          tmp1 <- unique(tmp1)
          for(row in 1:nrow(tmp1)){
            start <- tmp1[row,]$start-1
            cigarString <- tmp1[row,]$cigar
            count <- tmp1[row,]$count
            #Regular expression to split the CIGAR string into operations
            cigarSplit <- regmatches(cigarString, gregexpr("([0-9]+)([MIDNSHPX=])", cigarString))[[1]]
            numbers <- as.numeric(gsub("[A-Z=]", "", cigarSplit))
            characters <- gsub("[0-9]", "", cigarSplit)
            for(type in 1:length(characters)){
              if(characters[type] %in% c("M","X","=")){
                end=tail(start+1:numbers[type],n=1)
                positionCovered <- c(positionCovered,rep(paste0(start,":",end),count))
                start=end
              }else if(characters[type] %in% c("I")){
                start=start
              }else{
                start=tail(start+1:numbers[type],n=1)              }
            }
          }
          regionsCovered[[ref]] <- positionCovered
        }
        rm(tmp)
        
        #Filter covered bases for depth, calculate how many bases are covered per reference contig
        ##Function to give covered ranges
        expand_range <- function(range) {
          start_end <- unlist(strsplit(range, ":"))
          seq(as.numeric(start_end[1]), as.numeric(start_end[2]))
        }
        refCovered <- data.table()
        for(ref in names(regionsCovered)){
          print(ref)
          ranges <- regionsCovered[[ref]]
          ranges <- unlist(lapply(ranges, expand_range))
          #Remove all sites which are not covered at least x times by a read
          ranges <- data.table(table(ranges))
          ranges <- ranges[N>0] #depth limit in order to consider a reference nucleotide sufficiently covered for counting
          #Bases taxonomising reference nucleotides covered > x times
          bases <- sum(ranges$N)
          #Check if there long regions covered and separated by uncovered regions
          x <- sort(as.numeric(ranges$ranges))
          x <- split(x, cumsum(c(1, diff(x) != 1)))
          x <- sapply(x, function(seq) length(seq) >= 75) #specify region length in order to count a region as sufficiently covered
          x <- sum(x)
          #Count length of covered nucleotides
          refCovered <- rbind(refCovered,data.table("ntCovered"=length(unique(ranges$ranges)),
                                                    "basesTaxonomised"=bases,
                                                    "regionsCovered"=x,
                                                    "reference"=ref))
        }

        #Add genome size information, calculate covered nucleotides and total amount of matching bases per species (taxid)
        refCovered$taxid <- as.numeric(sapply(str_extract_all(refCovered$reference, "(?<=taxid\\|)\\d+"), unique))
        refCovered[, basesTaxonomised := sum(basesTaxonomised), by = taxid]
        refCovered[, ntCovered := sum(ntCovered), by = taxid]
        refCovered[, regionsCovered := sum(regionsCovered), by = taxid]
        refCovered <- unique(refCovered[,c("ntCovered","basesTaxonomised","regionsCovered","taxid")])
        
        #Merge information with sam table
        sam <- unique(sam[,c("taxid","nreads.sp","percentReads","all.nreads.taxonomised","all.bases.taxonomised")])
        sam <- merge(sam, refCovered, by="taxid",all.x=T)
        sam <- merge(sam, refSize[,c("taxid","size")], by="taxid")
        
        #Calculate % of read bases taxonomised
        sam$pBasesUnique <- 100*sam$ntCovered/sam$all.bases.taxonomised
        sam$pBasesTaxonomised <- 100*sam$basesTaxonomised/sam$all.bases.taxonomised
        
        #Calculate % of genome covered
        sam$pRefCovered <- 100*sam$ntCovered/sam$size
        
        #Assign mapq, read, position, pident, pmatch and %-abundance filtering columns
        sam[, 'pidentCut' := 95]
        sam[, 'pmatchCut' := 95]
        
        #Read diamond output file of the corresponding sample to extract the number of deduplicated reads for below normalisation
        file <- paste0("002_6_",sample,"_90.0BS-0.0aln_0.0idt_1e-10ev_NArm_0.1pMin-GE_cts.txt")
        dedup <- read.table(file.path(countsDir,file),header=F,sep='\t')$V15[1]
  
        #Assign number of deduplicated reads
        sam[, nreads.deduplicated := dedup]
  
        #Taxonomise taxids to avoid doing this in the filtering loop
        dfTax <- data.frame(getTaxonomy(unique(sam$taxid),'C:/Users/Michael/Desktop/ncbi_taxonomy_db/accessionTaxa.sql'))
        dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
    
        #Merge with sam
        sam <- sam[dfTax, on = "taxid"]
        
        #Associate samples with date and time based on sample name
        sam[, sample := sample]
        sam[dtSamplingAnalysis, date := i.date, on = "sample"]
    
        #Remove genera from labCtrl
        sam <- sam[!genus %in% c(labCtrl$Genus,"Homo")]
        
        #Export file
        saveRDS(sam,rdsName)
        write.csv(sam,csvName,row.names=F)
      }
    }
  }
}

#Evaluate the effect of pmatch as a threshold
dfTaxCalls <- data.frame()
for(f in grep(".filtered.rds",list.files(),value=T)){
  print(f)
  dfTaxCalls <- rbind(dfTaxCalls,
                      readRDS(f))
}

#Print histogram of covered regions
pdf("histo_regions-covered.pdf")
dfTaxCalls %>%
  dplyr::select(species,regionsCovered,sample) %>%
  unique() %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(regionsCovered=mean(regionsCovered)) %>%
  dplyr::ungroup() %>%
  unique() %>%
  pull(regionsCovered) %>%
  log10(.) %>%
  hist(.)
dev.off()

#Filter dfTaxCalls
dfTaxCalls %>% dplyr::group_by(sample) %>% dplyr::filter(percentReads>=0.05 & regionsCovered >= 10) %>% dplyr::ungroup() %>% dplyr::pull(species) %>% unique() %>% sort()
dfTaxCalls <- dfTaxCalls %>% dplyr::group_by(sample) %>% dplyr::filter(percentReads>=0.05 & regionsCovered >= 10) %>% dplyr::ungroup()

#Normalise to RPM per replicate
tmp <- dfTaxCalls %>%
  dplyr::group_by(sample,species) %>%
  dplyr::mutate(nreads.sp=sum(nreads.sp)) %>%
  dplyr::mutate(nreads.sp.norm=nreads.sp/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(!c(nreads.deduplicated)) %>%
  unique() %>%
  dplyr::select(!nreads.sp) %>%
  dplyr::select(sample,date,superkingdom,phylum,class,order,family,genus,species,nreads.sp.norm) %>% 
  unique()
write.csv(tmp,"replicate_date_RPM_per_species.csv")

#Normalise to RPM per date
dfTaxCalls <- dfTaxCalls %>%
  dplyr::group_by(date,species) %>%
  dplyr::mutate(nreads.sp=sum(nreads.sp)) %>%
  dplyr::mutate(nreads.sp.norm=nreads.sp/(sum(nreads.deduplicated)/1e6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(!c(sample,nreads.deduplicated)) %>%
  unique() %>%
  dplyr::select(!nreads.sp)

#Save dfTaxCalls
saveRDS(dfTaxCalls,"dfTaxCalls.rds")

########################################
#SHOW HOW PHIBASE SPECIES VARY OVER TIME: BWA MAPPING CUTOFF ANALYSOS
########################################
dir.create(file.path(analysisDir,"09_field-sampling-timecourse_phibase-species_bwa-analysis"))
setwd(file.path(analysisDir,"09_field-sampling-timecourse_phibase-species_bwa-analysis"))

#Read in taxonomy calls from BWA mapping and filter according to previously determined thresholds 
dfTaxCalls <- readRDS(file.path(analysisDir,"08_field-sampling-timecourse_phibase-species_bwa-cutoff-analysis","dfTaxCalls.rds"))

#Load phibase data for merging with dfTaxCalls (i.e. to add the host information to dfTaxCalls)
df <- read.csv(file.path(dataDir,"phi-base_v4-12_2021-09-02.csv")) %>% dplyr::select(Pathogen.ID,Host.ID)
colnames(df) <- c("taxid.pathogen","taxid.host")
#Call taxa on phibase data
df$taxid.pathogen <- as.numeric(df$taxid.pathogen) #needs to be numeric
dfTax <- data.frame(getTaxonomy(unique(df$taxid.pathogen),'C:/Users/Michael/Desktop/ncbi_taxonomy_db/accessionTaxa.sql'))
dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
#Control taxid content
setdiff(df$taxid.pathogen,dfTax$taxid)
setdiff(dfTax$taxid,df$taxid.pathogen)
length(intersect(dfTax$taxid,df$taxid.pathogen))
length(unique(df$taxid.pathogen))
#Merge to obtain phylogeny
df <- merge(df,dfTax,by.x="taxid.pathogen",by.y="taxid",all=T)
#Call taxa on phibase data
df$taxid.host <- as.numeric(df$taxid.host) #needs to be numeric
dfTax <- data.frame(getTaxonomy(unique(df$taxid.host),'C:/Users/Michael/Desktop/ncbi_taxonomy_db/accessionTaxa.sql'))
dfTax$taxid <- as.numeric(gsub(" ","",row.names(dfTax)))
#Control taxid content
setdiff(df$taxid.host,dfTax$taxid)
setdiff(dfTax$taxid,df$taxid.host)
length(intersect(dfTax$taxid,df$taxid.host))
length(unique(df$taxid.host))
#Merge to obtain phylogeny
df <- merge(df,dfTax,by.x="taxid.host",by.y="taxid",all=T)
#Exchange column names .x and .y
colnames(df) <- gsub(".x",".pathogen",colnames(df),fixed=T)
colnames(df) <- gsub(".y",".host",colnames(df),fixed=T)
#Filter taxonomised phibase dataframe for merging
df <- df %>% dplyr::select(superkingdom.pathogen,phylum.pathogen,class.pathogen,order.pathogen,family.pathogen,genus.pathogen,species.pathogen,superkingdom.host,phylum.host,class.host,order.host,family.host,genus.host,species.host)
#Add Puccinia triticina
df <- rbind(df,
            df %>% dplyr::select(superkingdom.pathogen,phylum.pathogen,class.pathogen,order.pathogen,family.pathogen,genus.pathogen,species.pathogen,superkingdom.host,phylum.host,class.host,order.host,family.host,genus.host,species.host) %>% dplyr::filter(species.pathogen=="Puccinia striiformis" & species.host=="Triticum aestivum") %>% unique() %>% dplyr::mutate(species.pathogen="Puccinia triticina"))
#Merge with hosts from phibase
colnames(dfTaxCalls) <- gsub("species","species.pathogen",colnames(dfTaxCalls))
dfTaxCalls <- dfTaxCalls %>% dplyr::filter(species.pathogen %in% df$species.pathogen) %>% dplyr::select(!c(superkingdom,phylum,class,order,family,genus))
dfTaxCalls <- merge(dfTaxCalls ,df,by="species.pathogen")

write.table(dfTaxCalls,"dfTaxCalls.phibasespecies.txt")

#Export table for supplementary
write.csv(dfTaxCalls %>% dplyr::select(nreads.sp.norm,species.pathogen,date) %>% unique() %>% reshape2::dcast(.,formula=species.pathogen~date,value.var="nreads.sp.norm",fun.aggregate=mean),"counts-over-time_SUPP.csv")

#Pathogen species found
length(setdiff(unique(dfTaxCalls$species.pathogen),NA))

#Host species targeted
length(setdiff(unique(dfTaxCalls$species.host),NA))

#Streptophyta pathogen species found
dfTaxCalls %>% dplyr::filter(phylum.host=="Streptophyta") %>% pull(species.pathogen) %>% unique() %>% length()

#Mammalia pathogen species found
dfTaxCalls %>% dplyr::filter(class.host=="Mammalia") %>% pull(species.pathogen) %>% unique() %>% length()

#H. sapiens pathogen species found
dfTaxCalls %>% dplyr::filter(species.host=="Homo sapiens") %>% pull(species.pathogen) %>% unique() %>% length()

#Show counts in treemaps and export counts as tables for all plant, non-plant and all pathogens
pdf("all_pathogens_per_pathogen-family_treemap.pdf")
treemap(dfTaxCalls %>%
          dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,phylum.host) %>%
          #dplyr::filter(phylum.host=="Streptophyta") %>%
          dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
          unique() %>%
          dplyr::group_by(species.pathogen) %>%
          dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::select(family.pathogen,cts,species.pathogen) %>%
          unique(), index = c("species.pathogen"), vSize = "cts", title = "", fontsize.labels = 25, fontface.labels = 4)
dev.off()
write.csv(dfTaxCalls %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,phylum.host) %>%
            #dplyr::filter(phylum.host=="Streptophyta") %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
            unique() %>%
            dplyr::group_by(species.pathogen) %>%
            dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(family.pathogen,cts,species.pathogen) %>%
            unique(),"all_pathogens_per_pathogen-family_treemap.csv")

pdf("plant_pathogens_per_pathogen-family_treemap.pdf")
treemap(dfTaxCalls %>%
          dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,phylum.host) %>%
          dplyr::filter(phylum.host=="Streptophyta") %>%
          dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
          unique() %>%
          dplyr::group_by(species.pathogen) %>%
          dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::select(family.pathogen,cts,species.pathogen) %>%
          unique(), index = c("species.pathogen"), vSize = "cts", title = "", fontsize.labels = 25, fontface.labels = 4)
dev.off()
write.csv(dfTaxCalls %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,phylum.host) %>%
            dplyr::filter(phylum.host=="Streptophyta") %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
            unique() %>%
            dplyr::group_by(species.pathogen) %>%
            dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(family.pathogen,cts,species.pathogen) %>%
            unique(),"plant_pathogens_per_pathogen-family_treemap.csv")

pdf("non-plant_pathogens_per_pathogen-family_treemap.pdf")
treemap(dfTaxCalls %>%
          dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,phylum.host) %>%
          dplyr::filter(!phylum.host=="Streptophyta") %>%
          dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
          unique() %>%
          dplyr::group_by(species.pathogen) %>%
          dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::select(family.pathogen,cts,species.pathogen) %>%
          unique(), index = c("species.pathogen"), vSize = "cts", title = "", fontsize.labels = 25, fontface.labels = 4)
dev.off()
write.csv(dfTaxCalls %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,phylum.host) %>%
            dplyr::filter(!phylum.host=="Streptophyta") %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
            unique() %>%
            dplyr::group_by(species.pathogen) %>%
            dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(family.pathogen,cts,species.pathogen) %>%
            unique(),"non-plant_pathogens_per_pathogen-family_treemap.csv")

#Only wheat, barley, pea pathogens
write.csv(dfTaxCalls %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm,species.host) %>%
            dplyr::filter(species.host %in% c("Triticum aestivum","Hordeum vulgare","Pisum sativum")) %>%
            dplyr::select(date,family.pathogen,species.pathogen,nreads.sp.norm) %>%
            unique() %>%
            dplyr::group_by(species.pathogen) %>%
            dplyr::mutate(cts=sum(nreads.sp.norm)) %>%
            dplyr::ungroup() %>%
            dplyr::select(family.pathogen,cts,species.pathogen) %>%
            unique(),"wheat_barley_pea-plant_pathogens_per_pathogen-family_treemap.csv")

#Create summary what pathogen species are colonising
write.csv(dfTaxCalls %>%
  dplyr::select(superkingdom.pathogen,phylum.pathogen,class.pathogen,order.pathogen,family.pathogen,genus.pathogen,species.pathogen,superkingdom.host,phylum.host,class.host,order.host,family.host,genus.host,species.host) %>%
  dplyr::group_by(phylum.host) %>%
  dplyr::mutate(SpeciesPerPhylum=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  unique(),"species.per.pathogen.csv")

dfTaxCalls %>%
  dplyr::select(species.pathogen,class.host) %>%
  unique() %>%
  dplyr::filter(class.host=="Mammalia") %>%
  pull(species.pathogen) %>%
  unique() %>%
  length()

dfTaxCalls %>%
  dplyr::select(species.pathogen,species.host) %>%
  unique() %>%
  dplyr::filter(species.host=="Homo sapiens") %>%
  pull(species.pathogen) %>%
  unique() %>%
  length()

#Show selected pathogen species by phylum and their abundance
pdf("bwamem_phibase-strptophyta-pathogens_treemap_phylum_species.pdf")
treemap(dfTaxCalls %>%
          dplyr::filter(phylum.host=="Streptophyta") %>%
          dplyr::select(date,species.pathogen,phylum.pathogen,nreads.sp.norm) %>%
          unique() %>%
          dplyr::group_by(date,species.pathogen,phylum.pathogen) %>%
          dplyr::mutate(nreads.sp.percent=100*nreads.sp.norm/sum(nreads.sp.norm)) %>%
          dplyr::ungroup() %>%
          dplyr::select(date,phylum.pathogen,species.pathogen,nreads.sp.percent) %>%
          unique() %>%
          dplyr::group_by(species.pathogen,phylum.pathogen) %>%
          dplyr::mutate(sumReads=sum(nreads.sp.percent)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sumReads=100*sumReads/sum(nreads.sp.percent)) %>%
          dplyr::select(phylum.pathogen,species.pathogen,sumReads) %>%
          unique(), index = c("phylum.pathogen","species.pathogen"), vSize = "sumReads", title = "", fontsize.labels = 25, fontface.labels = 4)
dev.off()

#Show which host has most pathogens
pdf("pathogens_per_host-family_treemap.pdf")
treemap(dfTaxCalls %>%
  dplyr::select(date,species.pathogen,family.host,phylum.host) %>%
  dplyr::filter(phylum.host=="Streptophyta") %>%
  unique() %>%
  dplyr::group_by(family.host) %>%
  dplyr::mutate(cts=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  dplyr::select(family.host,cts,species.pathogen) %>%
  unique(), index = c("family.host"), vSize = "cts", title = "", fontsize.labels = 25, fontface.labels = 4)
dev.off()

#Table of which host as most pathogens
write.csv(dfTaxCalls %>%
            dplyr::select(date,species.pathogen,family.host,phylum.host) %>%
            dplyr::filter(phylum.host=="Streptophyta") %>%
            unique() %>%
            dplyr::group_by(family.host) %>%
            dplyr::mutate(cts=length(unique(species.pathogen))) %>%
            dplyr::ungroup() %>%
            dplyr::select(family.host,cts,species.pathogen) %>%
            unique(),"pathogens_per_host_family.csv")

#Which plant families are affected by the detected pathogens
pdf("number-of-pathogens-per-host-family.pdf",h=5,w=12)
tmp <- df %>%
  dplyr::select(species.pathogen,family.host) %>%
  unique() %>% 
  dplyr::group_by(family.host) %>%
  dplyr::mutate(PHIbase=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  dplyr::select(family.host,PHIbase) %>%
  unique()%>%
  dplyr::arrange(PHIbase)

tmp <- dfTaxCalls %>% 
  dplyr::filter(phylum.host=="Streptophyta") %>%
  dplyr::select(species.pathogen,family.host) %>%
  unique() %>% 
  dplyr::group_by(family.host) %>%
  dplyr::mutate(Detected=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  dplyr::select(family.host,Detected) %>%
  unique() %>%
  merge(.,tmp,by="family.host") %>%
  dplyr::mutate(q=paste0(Detected,"/",PHIbase,"\n",round(100*Detected/PHIbase,0),"%")) %>%
  dplyr::filter(!family.host=="<NA>") %>%
  reshape2::melt(.)

order <- tmp %>% dplyr::filter(variable=="Detected") %>% dplyr::arrange(value) %>% pull(family.host) %>% rev()
ggplot() +
  geom_bar(data=tmp %>% dplyr::filter(variable=="PHIbase"),aes(fill=variable, y=value, x=factor(family.host,levels=order)), position='identity', stat="identity") +
  geom_bar(data=tmp %>% dplyr::filter(variable=="Detected"),aes(fill=variable, y=value, x=factor(family.host,levels=order)), position='identity', stat="identity") +
  geom_text(data=tmp %>% dplyr::filter(variable=="PHIbase"),aes(label=q, y=value, x=factor(family.host,levels=order)), vjust=-0.2, size=3) +
  theme_bw(base_size = 14) +
  scale_fill_aaas() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=16, face="italic"),
        strip.text = element_text(size = 10),
        legend.text=element_text(size=12),
        legend.title=element_text(size=16),
        legend.position="right") +
  labs(fill = "") +
  geom_hline(yintercept=85,alpha=-0) +
  xlab("") +
  ylab ("Pathogen species") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#How many pathogens do we find for the grown species nearby
pdf("number-of-pathogens-per-selected-host-species.pdf",h=5,w=4)
tmp <- df %>%
  dplyr::select(species.pathogen,species.host) %>%
  unique() %>% 
  dplyr::group_by(species.host) %>%
  dplyr::mutate(PHIbase=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  dplyr::select(species.host,PHIbase) %>%
  unique()%>%
  dplyr::arrange(PHIbase)

tmp <- dfTaxCalls %>% 
  dplyr::filter(species.host %in% c("Hordeum vulgare","Triticum aestivum","Pisum sativum")) %>%
  dplyr::select(species.pathogen,species.host) %>%
  unique() %>% 
  dplyr::group_by(species.host) %>%
  dplyr::mutate(Detected=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  dplyr::select(species.host,Detected) %>%
  unique() %>%
  merge(.,tmp,by="species.host") %>%
  dplyr::mutate(q=paste0(Detected,"/",PHIbase,"\n",round(100*Detected/PHIbase,0),"%")) %>%
  dplyr::filter(!species.host=="<NA>") %>%
  reshape2::melt(.)

order <- tmp %>% dplyr::filter(variable=="Detected") %>% dplyr::arrange(value) %>% pull(species.host) %>% rev()
ggplot() +
  geom_bar(data=tmp %>% dplyr::filter(variable=="PHIbase"),aes(fill=variable, y=value, x=factor(species.host,levels=order)), position='identity', stat="identity") +
  geom_bar(data=tmp %>% dplyr::filter(variable=="Detected"),aes(fill=variable, y=value, x=factor(species.host,levels=order)), position='identity', stat="identity") +
  geom_text(data=tmp %>% dplyr::filter(variable=="PHIbase"),aes(label=q, y=value, x=factor(species.host,levels=order)), vjust=-0.2, size=4) +
  theme_bw(base_size = 14) +
  scale_fill_aaas() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=16, face="italic"),
        strip.text = element_text(size = 10),
        legend.text=element_text(size=12),
        legend.title=element_text(size=16),
        legend.position="right") +
  labs(fill = "") +
  geom_hline(yintercept=25,alpha=-0) +
  xlab("") +
  ylab ("Pathogen species") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()

#Count which pathogens are particularly abundant by host and make read plots over time
pdf("sumReads-of-pathogens.pdf",h=10,w=8)
dfTaxCalls %>%
  dplyr::select(date,species.pathogen,nreads.sp.norm,species.host,family.host,phylum.host) %>%
  unique() %>%
  dplyr::filter(phylum.host=="Streptophyta") %>%
  dplyr::select(family.host,species.pathogen,family.host,date,nreads.sp.norm,date) %>%
  unique() %>%
  dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>%
  ggplot() +
  geom_tile(aes(x=as.character(dateReformat), y=factor(species.pathogen,levels=sort(unique(species.pathogen),decreasing=T)), fill=nreads.sp.norm)) +
  coord_equal() +
  scale_fill_gradientn(trans='log10', colours = rev(brewer.pal(11, 'RdYlBu'))) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=6),
        axis.text=element_text(size=6),
        axis.title=element_text(size=6),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=5),
        axis.text.y = element_text(face = "italic"),
        legend.text=element_text(size=13),
        legend.title=element_text(size=12),
        legend.position="right") +
  labs(fill = "RPM dedup.") +
  #facet_wrap(~species.host,scales="free_y") +
  xlab("") +
  ylab ("") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()


#Plot number of species per family that are colonised per pathogen
tmp <- dfTaxCalls %>%
  dplyr::filter(phylum.host=="Streptophyta") %>%
  dplyr::select(family.host,species.pathogen,species.host) %>%
  dplyr::group_by(family.host,species.pathogen) %>%
  dplyr::mutate(count=length(unique(species.host))) %>%
  dplyr::ungroup() %>%
  dcast(.,formula=species.pathogen~family.host,value.var="count",fun.aggregate=mean) %>%
  reshape2::melt(.,id.var="species.pathogen") %>%
  dplyr::filter(!is.na(value)) %>%
  unique()

#Rbind tmp with sum of pathogen species
tmp <- rbind(tmp, tmp %>%
               dplyr::group_by(species.pathogen) %>%
               dplyr::mutate(value=sum(value)) %>%
               dplyr::ungroup() %>%
               dplyr::mutate(variable="SUM") %>%
               unique())
#Define order of hosts based on number of pathogens that colonise them
order.hosts <- tmp %>% group_by(variable) %>%
  dplyr::mutate(value=length(unique(species.pathogen))) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable,value) %>%
  unique() %>%
  dplyr::arrange(desc(value)) %>%
  pull(variable)
#Custom color palette for plotting
colorpalette <- c(brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'),brewer.pal(12,'Paired'),brewer.pal(9,'Set1'))
#Plot
pdf("host_families_number_species_per-pathogen.pdf",width=5)
tmp %>%
  ggplot() +
  geom_point(aes(x=factor(variable,levels=order.hosts),
                 y=factor(species.pathogen,levels=sort(unique(species.pathogen),decreasing=T)),
                 size=log10(value), colour=variable),shape=21) +
  theme_bw(base_size = 14) +
  scale_color_manual(values=colorpalette) +
  theme(plot.title = element_text(hjust = 0.5, size=6),
        axis.text=element_text(size=6),
        axis.title=element_text(size=6),
        axis.text.x=element_text(angle=60, hjust = 1.0,size=6),
        axis.text.y = element_text(face = "italic"),
        legend.text=element_text(size=13),
        legend.title=element_text(size=12),
        legend.position="") +
  labs(fill = "") +
  xlab("") +
  ylab ("") +
  ggtitle("") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off()


pdf("sumReads-of-pathogens_wheat_barley_pea.pdf",h=12,w=10)
#Calculate the lower and upper end of the scaling factor for the tile plots
max.scaling <- 0
for(host in c("Triticum aestivum","Hordeum vulgare","Pisum sativum")){
  n <- dfTaxCalls %>% dplyr::filter(species.host==host) %>% pull(nreads.sp.norm) %>% max() %>% ceiling()
  ifelse(n>max.scaling,max.scaling<-n,"")
}
min.scaling <- 1e12
for(host in c("Triticum aestivum","Hordeum vulgare","Pisum sativum")){
  n <- dfTaxCalls %>% dplyr::filter(species.host==host) %>% pull(nreads.sp.norm) %>% min() %>% floor()
  ifelse(n<min.scaling,min.scaling<-n,"")
}

#Define good breaks for the colour legend
print(min.scaling)
print(max.scaling)
plot.breaks <- c(10,100,1000,10000)

for(host in c("Triticum aestivum","Hordeum vulgare","Pisum sativum")){
  allDatesReformat <- dfTaxCalls %>% dplyr::arrange(date) %>% dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>% pull(dateReformat) %>% unique()
  tmp <- dfTaxCalls %>%
    dplyr::filter(species.host==host) %>%
    dplyr::select(date,species.pathogen,nreads.sp.norm,species.host,family.host,phylum.host) %>%
    unique() %>%
    dplyr::select(species.host,species.pathogen,date,nreads.sp.norm,date) %>%
    unique()
  order <- tmp %>% dplyr::group_by(species.pathogen) %>% dplyr::mutate(m=sum(nreads.sp.norm)) %>% dplyr::ungroup() %>% dplyr::select(species.pathogen,m) %>% unique() %>% dplyr::arrange(m) %>% dplyr::pull(species.pathogen)
  tmp %>%
    dplyr::mutate(dateReformat=gsub(".2015","",format(date,format="%d.%m.%Y"))) %>%
    ggplot() +
    geom_tile(aes(x=as.character(dateReformat), y=factor(species.pathogen,levels=order), fill=nreads.sp.norm)) +
    coord_equal() +
    scale_fill_gradientn(trans='log10', colours = rev(brewer.pal(11, 'RdYlBu')), limits=c(min.scaling,max.scaling), breaks=plot.breaks) +
    theme(plot.title = element_text(hjust = 0.5, size=6),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          axis.text.x=element_text(angle=60, hjust = 1.0,size=14),
          axis.text.y = element_text(face = "italic"),
          legend.text=element_text(size=15),
          legend.title=element_text(size=12),
          strip.text = element_text(size = 16, face = "italic"),
          legend.position="right") +
    labs(fill = "RPM dedup.") +
    facet_wrap(~species.host) +
    xlab("") +
    ylab ("") +
    ggtitle("") +
    scale_x_discrete(limits=allDatesReformat) +
    guides(color=guide_legend(title="")) -> p
  plot(p)
}
dev.off()

########################################
#PLOT WEATHER VARIABLES
########################################
dir.create(file.path(analysisDir,"10_field-sampling-timecourse_weather"))
setwd(file.path(analysisDir,"10_field-sampling-timecourse_weather"))

##Filter date table by min/max
dfWeather <- dfWeather %>% dplyr::filter(date %in% as.Date(seq(min(dfTaxCalls$date), max(dfTaxCalls$date),by="day")))

#Extract maximal weather values and plot to associate with species
for(measure in setdiff(colnames(dfWeather)[seq(4,19)],"Wind.Direction")){
  print(measure)
  pdf(paste0("MAX_weather_variable_",measure,".pdf"))
  tmp <- dfWeather %>%
    dplyr::select(measure,date)
  colnames(tmp) <- c("measure","date")
  tmp <- tmp %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(measure=max(measure)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date %in% seq(as.Date("2015-06-08"),as.Date("2015-07-29"),by="day")) %>%
    unique()
  tmp <- rbind(tmp,data.frame(date=seq(as.Date("2015-06-08"),as.Date("2015-06-16"),by="day"),"measure"=NaN))
  tmp %>%
    ggplot() +
    geom_tile(aes(x=as.Date(date), y="1", fill=measure)) +
    coord_equal() +
    scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdYlBu'))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.text=element_text(size=10),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=0, hjust = 1.0,size=20),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          legend.position="right") +
    labs(fill = "X") +
    #facet_wrap(~species.host,scales="free_y") +
    xlab("") +
    ylab ("") +
    ggtitle(measure) +
    guides(color=guide_legend(title="")) -> p
  plot(p)
  tmp %>%
    ggplot() +
    geom_tile(aes(x=as.character(date), y="1", fill=measure)) +
    coord_equal() +
    scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdYlBu'))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.text=element_text(size=10),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=60, hjust = 1.0,size=8),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          legend.position="right") +
    labs(fill = "X") +
    #facet_wrap(~species.host,scales="free_y") +
    xlab("") +
    ylab ("") +
    ggtitle(measure) +
    guides(color=guide_legend(title="")) -> p
  plot(p)
  dev.off()
}

#Extract mean weather values and plot to associate with species
for(measure in setdiff(colnames(dfWeather)[seq(4,19)],"Wind.Direction")){
  print(measure)
  pdf(paste0("MEAN_weather_variable_",measure,".pdf"))
  tmp <- dfWeather %>%
    dplyr::select(measure,date)
  colnames(tmp) <- c("measure","date")
  tmp <- tmp %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(measure=mean(measure)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date %in% seq(as.Date("2015-06-08"),as.Date("2015-07-29"),by="day")) %>%
    unique()
  tmp <- rbind(tmp,data.frame(date=seq(as.Date("2015-06-08"),as.Date("2015-06-16"),by="day"),"measure"=NaN))
  tmp %>%
    ggplot() +
    geom_tile(aes(x=as.Date(date), y="1", fill=measure)) +
    coord_equal() +
    scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdYlBu'))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.text=element_text(size=10),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=0, hjust = 1.0,size=20),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          legend.position="right") +
    labs(fill = "X") +
    #facet_wrap(~species.host,scales="free_y") +
    xlab("") +
    ylab ("") +
    ggtitle(measure) +
    guides(color=guide_legend(title="")) -> p
  plot(p)
  tmp %>%
    ggplot() +
    geom_tile(aes(x=as.character(date), y="1", fill=measure)) +
    coord_equal() +
    scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdYlBu'))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.text=element_text(size=10),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=60, hjust = 1.0,size=8),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          legend.position="right") +
    labs(fill = "X") +
    #facet_wrap(~species.host,scales="free_y") +
    xlab("") +
    ylab ("") +
    ggtitle(measure) +
    guides(color=guide_legend(title="")) -> p
  plot(p)
  dev.off()
}

#Plot wind direction
pdf("wind_direction.pdf")
tmp <- dfSampling %>%
  dplyr::filter(location == "field" & collection.minutes == 60) %>%
  dplyr::mutate(date=as.Date(date, format = "%d/%m/%Y")) %>%
  dplyr::select(wind.direction,date) %>%
  dplyr::filter(date %in% seq(as.Date("2015-06-08"),as.Date("2015-07-29"),by="day")) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(wind.direction=ifelse(length(unique(wind.direction))>1,"Variable",wind.direction)) %>%
  dplyr::ungroup() %>%
  unique()
colnames(tmp) <- c("measure","date")
#Order based on direction
tmp$measure <- factor(tmp$measure, levels=c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "Variable", "No wind"))

tmp %>%
  ggplot() +
  geom_tile(aes(x=as.Date(date), y="1", fill=measure)) +
  coord_equal() +
  scale_fill_aaas() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=0, hjust = 1.0,size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.position="right") +
  labs(fill = "X") +
  #facet_wrap(~species.host,scales="free_y") +
  xlab("") +
  ylab ("") +
  ggtitle("Wind dir.") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off() 

#Plot sampling weather
pdf("sampling_weather.pdf")
tmp <- dfSampling %>%
  dplyr::filter(location == "field" & collection.minutes == 60) %>%
  dplyr::mutate(date=as.Date(date, format = "%d/%m/%Y")) %>%
  dplyr::select(sampling.weather,date) %>%
  dplyr::filter(date %in% seq(as.Date("2015-06-08"),as.Date("2015-07-29"),by="day")) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(sampling.weather=ifelse(length(unique(sampling.weather))>1,"Variable",sampling.weather)) %>%
  dplyr::ungroup() %>%
  unique()
colnames(tmp) <- c("measure","date")
#Order based on direction
tmp$measure <- factor(tmp$measure, levels=c("Rain", "Dry", "Mist",  "Variable"))

tmp %>%
  ggplot() +
  geom_tile(aes(x=as.Date(date), y="1", fill=measure)) +
  coord_equal() +
  scale_fill_aaas() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=0, hjust = 1.0,size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.position="right") +
  labs(fill = "X") +
  #facet_wrap(~species.host,scales="free_y") +
  xlab("") +
  ylab ("") +
  ggtitle("Sampling weather") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off() 


#Plot plot from which wind came
pdf("wind_direction_plot.pdf")
tmp <- dfSampling %>%
  dplyr::filter(location == "field" & collection.minutes == 60) %>%
  dplyr::mutate(date=as.Date(date, format = "%d/%m/%Y")) %>%
  dplyr::select(wind.plot,date) %>%
  dplyr::filter(date %in% seq(as.Date("2015-06-08"),as.Date("2015-07-29"),by="day")) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(wind.plot=ifelse(length(unique(wind.plot))>1,"Variable",wind.plot)) %>%
  dplyr::ungroup() %>%
  unique()
colnames(tmp) <- c("measure","date")
#Order based on direction
tmp$measure <- factor(tmp$measure, levels=c("Wheat", "Wheat, Barley", "Barley", "Pea", "Variable", "No wind"))

tmp %>%
  ggplot() +
  geom_tile(aes(x=as.Date(date), y="1", fill=measure)) +
  coord_equal() +
  scale_fill_aaas() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=0, hjust = 1.0,size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.position="right") +
  labs(fill = "X") +
  #facet_wrap(~species.host,scales="free_y") +
  xlab("") +
  ylab ("") +
  ggtitle("Plot direction") +
  guides(color=guide_legend(title="")) -> p
plot(p)
dev.off() 

#Match field-samples with wind speed
tmp <- dfSampling %>%
  dplyr::filter(location == "field" & collection.minutes == 60) %>%
  dplyr::select(sample,date,time.start)
tmp1 <- data.frame()
for(s in tmp$sample){
  day <- tmp %>% dplyr::filter(sample==s) %>% pull(date) %>% as.Date(., format = "%d/%m/%Y")
  time.start <- tmp %>% dplyr::filter(sample==s) %>% pull(time.start) %>% as.ITime(.)
  if(!is.na(time.start) & day %in% dfWeather$date){
    time.end <- as.ITime(time.start) + as.ITime("01:00:00")
    speed <- dfWeather %>% filter(date==day & as.ITime(Time)>=time.start & as.ITime(Time)<=time.end) %>% pull(Wind.Speed.km.h.) %>% mean()
    tmp1 <- rbind(tmp1,
                  data.frame("sample"=s,"day"=day,"avg.speed"=speed))
  }else{
    tmp1 <- rbind(tmp1,
                  data.frame("sample"=s,"day"=day,"avg.speed"="Not available"))
  }
}
write.csv(tmp1,"average_wind-speed_sampling.csv")

########################################
#GROUP AND PLOT PLOT PATHOGENS AS LINEPLOTS VARIABLES
########################################
dir.create(file.path(analysisDir,"11_field-sampling-timecourse_lineplots_correlation"))
setwd(file.path(analysisDir,"11_field-sampling-timecourse_lineplots_correlation"))

#Correlate normalised species counts with weather data
dfWeatherPathogen <- merge(dfWeather %>% dplyr::mutate(date=as.Date(Time)) %>% dplyr::select(!c("NO.","Time","Interval","Wind.Direction")) %>% reshape2::melt(.,id.vars="date") %>% dplyr::group_by(date,variable) %>% dplyr::mutate(value=mean(value)) %>% dplyr::ungroup() %>% unique(),
                           dfTaxCalls %>% dplyr::filter(species.host %in% c("Triticum aestivum","Hordeum vulgare","Pisum sativum")),
                           by="date")

dfCorrelation <- data.frame()
for(species in unique(dfWeatherPathogen$species.pathogen)){
  tmp <- dfWeatherPathogen %>% dplyr::filter(species.pathogen==species) %>% dplyr::select(variable,value,date,nreads.sp.norm,species.pathogen) %>% unique()
  if(length(unique(tmp$date))>5){
    dfCorrelation <- rbind(dfCorrelation,
                           tmp %>% dplyr::arrange(date) %>% dplyr::group_by(variable) %>% dplyr::mutate(R=cor(nreads.sp.norm,value,method="pearson")) %>% dplyr::ungroup())
  }
}

#Select species which correlate with at least one variable well
df <- dfCorrelation %>%
  dplyr::filter(variable %in% c("Outdoor.Temperature.oC.","Outdoor.Humidity...","Wind.Speed.km.h.","X24.Hour.Rainfall.mm.")) %>%
  dplyr::mutate(variable=gsub("Outdoor.Temperature.oC.","Temperature [C]",variable)) %>%
  dplyr::mutate(variable=gsub("Outdoor.Humidity...","Humidity [%]",variable)) %>%
  dplyr::mutate(variable=gsub("Wind.Speed.km.h.","Wind Speed [km/h]",variable)) %>%
  dplyr::mutate(variable=gsub("X24.Hour.Rainfall.mm.","Rainfall 24h [mm]",variable)) %>%
  dplyr::select(variable,R,species.pathogen) %>%
  unique() %>%
  reshape2::dcast(.,variable~species.pathogen,value.var="R")
row.names(df) <- df$variable
df <- df %>% dplyr::select(!variable) 
df <- as.matrix(t(df))

#Heatmap of correlation: Weather variables with normalised counts
x <- pheatmap(df,display_numbers=T)
pdf("heatmap_species_weather.pdf",height=4,width=5)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()

#Cluster patterns of pathogens
set.seed(123)
tmp <- dfTaxCalls %>%
  dplyr::filter(phylum.host %in% c("Streptophyta")) %>%
  dplyr::select(date,species.pathogen,nreads.sp.norm) %>%
  unique() %>%
  dplyr::group_by(species.pathogen) %>%
  dplyr::mutate(scaled=scales::rescale(nreads.sp.norm)) %>%
  dplyr::ungroup() %>%
  dplyr::select(date,species.pathogen,scaled) %>%
  unique()
tmp <- reshape2::dcast(species.pathogen~date,data=tmp)
tmp[is.na(tmp)] <- 0
row.names(tmp) <- tmp$species.pathogen
tmp <- tmp %>% dplyr::select(!species.pathogen)
dfAP <- apcluster(s=corSimMat(tmp, sel=NA, r=1, signed=TRUE, method="pearson"))
#Extract a data frame containing DE gene and cluster information
dfClusters <- data.frame()
for(cluster in seq(1,length(dfAP@clusters))){
  dfClusters <- rbind(dfClusters, data.frame("species.pathogen"=names(dfAP[[cluster]]), "cluster"=cluster))
}
write.csv(dfClusters %>% merge(.,dfTaxCalls %>% dplyr::select(!c(date,nreads.sp.norm)) %>% unique(),by="species.pathogen"),"dfClusters.csv")

#Plot plant pathogen data alongside with weather fit 
pdf("species_plots_geomline_all-clusters.pdf",h=10,w=15)
tmp <- dfTaxCalls %>%
  dplyr::filter(phylum.host=="Streptophyta") %>%
  dplyr::select(date,species.pathogen,nreads.sp.norm,species.host) %>%
  unique() %>%
  dplyr::group_by(species.pathogen) %>%
  dplyr::mutate(scaled=scales::rescale(nreads.sp.norm)) %>%
  dplyr::ungroup() %>%
  unique() %>% 
  merge(.,dfClusters,by="species.pathogen") %>% 
  dplyr::group_by(cluster,date) %>%
  dplyr::mutate(mean.scaled=mean(scaled)) %>%
  dplyr::ungroup()

ggplot() +
  geom_point(data=tmp, aes(x=date, y=scaled, group=interaction(species.pathogen,date),colour=species.host),shape=21,colour="black",size=0.5,lwd=0.5) +
  geom_line(data=tmp, aes(x=date, y=scaled, group=species.pathogen,colour=species.host),alpha=0.5) +
  geom_line(data=tmp %>% dplyr::select(date,cluster,mean.scaled) %>% unique(), aes(x=date, y=mean.scaled),linetype="dashed",colour="grey34",size=1) +
  #scale_fill_gradientn(trans='log10',colours = rev(brewer.pal(11, 'RdYlBu'))) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=90, size=20),
        strip.text = element_text(size = 15),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position="") +
  labs(colour = "Host") +
  facet_wrap(~cluster) +
  xlab("") +
  ylab("Scaled, RPM dedup. [0-1]") +
  ggtitle("") +
  guides(alpha='none') -> p
plot(p)
dev.off()

pdf("species_plots_geomline_all-clusters_selected-species-coloured.pdf",h=5,w=10)
tmp <- dfTaxCalls %>%
  dplyr::filter(phylum.host=="Streptophyta") %>%
  dplyr::select(date,species.pathogen,nreads.sp.norm,species.host) %>%
  unique() %>%
  dplyr::group_by(species.pathogen) %>%
  dplyr::mutate(scaled=scales::rescale(nreads.sp.norm)) %>%
  dplyr::ungroup() %>%
  unique() %>% 
  merge(.,dfClusters,by="species.pathogen") %>% 
  dplyr::group_by(cluster,date) %>%
  dplyr::mutate(mean.scaled=mean(scaled)) %>%
  dplyr::ungroup() %>%
  #Assign colouring variable by host
  dplyr::group_by(species.host) %>%
  dplyr::mutate(host.colouring=ifelse(species.host %in% c("Triticum aestivum","Hordeum vulgare","Pisum sativum"),species.host,"Other")) %>%
  dplyr::mutate(host.alpha=ifelse(species.host %in% c("Triticum aestivum","Hordeum vulgare","Pisum sativum"),1,0.90)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(cluster %in% c(1,6)) %>%
  dplyr::select(!species.host) %>%
  unique()

ggplot() +
  #geom_point(aes(x=date, y=scaled, group=interaction(species.pathogen,date),fill=nreads.sp.norm)shape=21,colour="black",size=3) +
  #geom_line(data=tmp,aes(x=date, y=scaled, group=species.pathogen,colour=host.colouring),alpha=0.5) +
  geom_line(data=tmp %>% dplyr::filter(host.colouring=="Other"), aes(x=date, y=scaled, group=species.pathogen,colour=factor(host.colouring,levels=c("Triticum aestivum","Hordeum vulgare","Pisum sativum","Other")), alpha=as.factor(host.alpha))) +
  geom_line(data=tmp %>% dplyr::filter(!host.colouring=="Other"), aes(x=date, y=scaled, group=species.pathogen,colour=factor(host.colouring,levels=c("Triticum aestivum","Hordeum vulgare","Pisum sativum","Other"))),alpha=0.6) +
  #geom_line(data=tmp %>% dplyr::select(date,cluster,mean.scaled) %>% unique(), aes(x=date, y=mean.scaled),linetype="dashed",colour="grey34",size=1) +
  geom_smooth(data=tmp, aes(x=date, y=scaled),linetype="dashed",colour="grey34",size=1,se=T) +
  #scale_fill_gradientn(trans='log10',colours = rev(brewer.pal(11, 'RdYlBu'))) +
  scale_color_aaas() +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust=1, size=20),
        strip.text = element_text(size = 15),
        legend.text=element_text(size=18, face = "italic"),
        legend.title=element_text(size=20),
        legend.position="right") +
  labs(colour = "Host") +
  facet_wrap(~cluster,nrow=1) +
  xlab("") +
  ylab("Scaled, RPM dedup. [0-1]") +
  ggtitle("") +
  scale_y_continuous(breaks=c(0,0.5,1.0)) +
  guides(alpha='none') -> p
plot(p)
dev.off()

#Plot each species as geom_line plot
pdf("species_plots_geomline_selected-species.pdf",h=5,w=5)
for(sp in dfTaxCalls %>% dplyr::filter(species.host %in% c("Triticum aestivum","Hordeum vulgare")) %>% dplyr::pull(species.pathogen) %>% unique()){
  tmp <- dfTaxCalls %>%
    dplyr::select(date,species.pathogen,nreads.sp.norm,species.host) %>%
    unique() %>%
    unique() %>%
    dplyr::group_by(species.pathogen) %>%
    dplyr::mutate(scaled=scales::rescale(nreads.sp.norm)) %>%
    dplyr::ungroup() %>%
    unique() %>% 
    merge(.,dfClusters,by="species.pathogen") %>% 
    dplyr::group_by(cluster,date) %>%
    dplyr::mutate(mean.scaled=mean(scaled)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(species.pathogen==sp)
  ggplot() +
    geom_point(data=tmp, aes(x=date, y=scaled, group=interaction(species.pathogen,date),fill=nreads.sp.norm),shape=21,colour="black",size=3) +
    geom_line(data=tmp, aes(x=date, y=scaled, group=species.pathogen),colour="grey34",linetype="dashed") +
    geom_smooth(data=tmp, aes(x=date, y=scaled, group=species.pathogen),colour="grey34",linetype="dashed") +
    scale_fill_gradientn(trans='log10',colours = rev(brewer.pal(11, 'RdYlBu'))) +
    scale_color_npg() +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=20, face = "italic"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=60, hjust=1, size=20),
          strip.text = element_text(size = 15),
          legend.text=element_text(size=18),
          legend.title=element_text(size=20),
          legend.position="right") +
    labs(fill = "RPM dedup.") +
    scale_y_continuous(breaks=c(0.0,0.5,1.0)) +
    xlab("") +
    ylab("Scaled, RPM dedup. [0-1]") +
    ggtitle(sp) +
    guides(alpha='none') -> p
  plot(p)
}
dev.off()

########################################
#FIND SNPS
########################################
dir.create(file.path(analysisDir,"12_PST-SNP-detection"))
setwd(file.path(analysisDir,"12_PST-SNP-detection"))

#Read in taxonomy calls from BWA mapping and filter according to previously determined thresholds 
dfTaxCalls <- readRDS(file.path(analysisDir,"08_field-sampling-timecourse_phibase-species_bwa-cutoff-analysis","dfTaxCalls.rds"))

#Load run accessions and sample id data
fieldIsolates <- c("13/27","13/21","13/14","13/22","13/65","13/12","13/09","13/29","13/15","13/25","13/34","13/19","13/26","13/23","13/24","13/32","13/20","13/30","13/38","13/28","13/33","13/36","13/39","13/71","13/42","13/40","13/35","CL1","13/18","13/37","13/520","13/123","T13/1","T13/2","T13/3","13/120","13/182","RB1","RB2")
dfPRJNA256347 <- read.table(file.path(dataDir,"filereport_read_run_PRJNA256347.tsv"),header=T,sep='\t') %>%
  dplyr::filter(sample_alias %in% fieldIsolates)

#Read reference tstv ratios for various depth and quality settings
dfTSTVref <- data.frame()
files <- grep("SRR",grep("filtered.calls.TSTV.txt",list.files(file.path(dataDir,"bwa_airseq_vs_PST")),value=T),value=T)
for(f in files){
  print(f)
  x <- read.table(file.path(dataDir,"bwa_airseq_vs_PST",f))
  if(nrow(x)>0){
    dfTSTVref <- rbind(dfTSTVref,
                          data.frame("sample"=strsplit(f,".",fixed=T)[[1]][1],
                                     "tstv"=x$V5,
                                     "qual"=unlist(strsplit(f,".",fixed=T))[3],
                                     "depth"=unlist(strsplit(f,".",fixed=T))[5]))
  }
}

dfTSTVref <- dfTSTVref %>%
  dplyr::filter(sample %in% dfPRJNA256347$run_accession) %>%
  dplyr::mutate(meanTSTV=mean(tstv)) %>%
  dplyr::mutate(sdTSTV=sd(tstv)) %>%
  dplyr::select(meanTSTV,sdTSTV) %>%
  unique() %>%
  dplyr::mutate("dataset"="reference")

#Dates when PST was found
datesPST <- dfTaxCalls %>% dplyr::filter(species=="Puccinia striiformis") %>% dplyr::pull(date) %>% unique() %>% as.Date()

#Read airseq tstv ratios
dfTSTVairseq <- data.frame()
files <- grep("airseq_JIC",grep("filtered.calls.TSTV.txt",list.files(file.path(dataDir,"bwa_airseq_vs_PST")),value=T),value=T)
for(f in files){
  print(f)
  x <- read.table(file.path(dataDir,"bwa_airseq_vs_PST",f))
  if(nrow(x)>0){
    dfTSTVairseq <- rbind(dfTSTVairseq,
                          data.frame("sample"=strsplit(f,".",fixed=T)[[1]][1],
                                     "tstv"=x$V5))
  }
}
dfTSTVairseq <- merge(dfTSTVairseq,dfSamplingAnalysis[c("sample","date")],by="sample")
dfTSTVairseq$date <- as.Date(dfTSTVairseq$date, format = "%d/%m/%Y")

dfTSTVairseq <- dfTSTVairseq %>%
  dplyr::filter(date %in% datesPST) %>%
  dplyr::mutate(meanTSTV=mean(tstv)) %>%
  dplyr::mutate(sdTSTV=sd(tstv)) %>%
  dplyr::select(meanTSTV,sdTSTV) %>%
  unique() %>%
  dplyr::mutate("dataset"="airseq")

write.csv(rbind(dfTSTVairseq,dfTSTVref),"tstv.csv")

#Read airseq SNP calls
dfSNPairseq <- data.frame()
files <- grep("airseq_JIC",grep("filtered.calls.vcf",list.files(file.path(dataDir,"bwa_airseq_vs_PST")),value=T),value=T)
if(!file.exists("dfSNPairseq.rds")){
  for(f in files){
  print(f)
  x <- readLines(file.path(dataDir,"bwa_airseq_vs_PST",f))
  if(!length(grep("#",x)) == length(x)){
    tmp <- read.table(file.path(dataDir,"bwa_airseq_vs_PST",f)) %>%
      dplyr::mutate(location=paste0(V1,"_",V2)) %>%
      dplyr::select(location,nt=V5) %>%
      dplyr::mutate(sample=strsplit(f,".",fixed=T)[[1]][1]) %>%
      dplyr::filter(nchar(nt)==1)
    dfSNPairseq <- rbind(dfSNPairseq,tmp)
    }
  }
  saveRDS(dfSNPairseq,"dfSNPairseq.rds")
}
dfSNPairseq <- readRDS("dfSNPairseq.rds")

#Filter the airseq data for the date when PST was called with all cutoffs
dfSNPairseq <- dfSNPairseq %>%
  merge(.,dfSamplingAnalysis[c("sample","date")],by="sample") %>%
  dplyr::mutate(date=as.Date(date, format = "%d/%m/%Y")) %>%
  dplyr::filter(date %in% datesPST)

#Loci and changes
length(unique(paste0(dfSNPairseq$location)))
length(unique(paste0(dfSNPairseq$location,dfSNPairseq$nt)))

#Control for overlaps with the loci of cluster 1-4 described by Hubbard et al.
lociOfInterest <- read.csv(file.path(dataDir,"supplementary_locus_table_provided.csv"),header=T,fileEncoding='UTF-8-BOM') %>%
  dplyr::mutate(location=paste0(contig,"_",position)) %>%
  dplyr::select(location,aa.cluster1=consensus_aa,aa.cluster2=consensus_aa.1,aa.cluster3=consensus_aa.2,aa.cluster4=consensus_aa.3)
tmp <- merge(dfSNPairseq,
             lociOfInterest,
             by="location")
dim(tmp)

#Read reference SNP calls
if(!file.exists("dfSNPreference.rds")){
  dfSNPreference <- data.frame()
  files <- grep("SRR",grep("filtered.calls.vcf",list.files(file.path(dataDir,"bwa_airseq_vs_PST")),value=T),value=T)
  for(f in files){
    print(f)
    x <- readLines(file.path(dataDir,"bwa_airseq_vs_PST",f))
    if(!length(grep("#",x)) == length(x)){
      tmp <- read.table(file.path(dataDir,"bwa_airseq_vs_PST",f)) %>%
        dplyr::mutate(location=paste0(V1,"_",V2)) %>%
        dplyr::select(location,nt=V5) %>%
        dplyr::mutate(n.loci=length(unique(paste0(location,nt)))) %>%
        dplyr::mutate(sample=strsplit(f,".",fixed=T)[[1]][1]) %>%
        merge(.,dfPRJNA256347 %>% dplyr::select(sample=run_accession,sample_alias),by="sample") %>%
        dplyr::filter(nchar(nt)==1)
      dfSNPreference <- rbind(dfSNPreference, tmp)
    }
  }
  dfSNPreference <- dfSNPreference %>% dplyr::filter(sample_alias %in% fieldIsolates)
  saveRDS(dfSNPreference,"dfSNPreference.rds")
}
dfSNPreference <- readRDS("dfSNPreference.rds")

#Select exclusively unique SNPs per isolate for the reference data
dfSNPreference <- dfSNPreference %>%
  dplyr::group_by(location,nt) %>%
  dplyr::filter(length(unique(sample))==1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n.loci.unique.reference=length(unique(paste0(location,nt)))) %>%
  dplyr::ungroup() %>%
  unique()
write.csv(dfSNPreference,"unique_snps.csv")
length(unique(paste0(dfSNPreference$location,dfSNPreference$nt)))

#Which strain is most present in the air-seq data?
tmp <- dfSNPairseq %>%
  dplyr::mutate(n.loci.airseq=length(unique(paste0(location,nt)))) %>%
  dplyr::select(location,nt,n.loci.airseq,sample.airseq=sample) %>%
  unique() %>%
  merge(.,
        dfSNPreference,
        by=c("location","nt")) %>%
  dplyr::group_by(sample_alias) %>%
  dplyr::mutate(n.loci.shared=length(unique(location))) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample,n.loci.shared,n.loci.airseq,sample_alias,n.loci.unique.reference) %>%
  unique()
write.csv(tmp,"SNP_overlap.summarised.csv")

#Control the overlap between the time-filtered airseq data and the reference data filtered for unique SNPs
tmp <- dfSNPairseq %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(n.loci.airseq.per.date=length(unique(paste0(location,nt)))) %>%
  dplyr::ungroup() %>%
  dplyr::select(location,nt,date,n.loci.airseq.per.date,sample.airseq=sample) %>%
  unique() %>%
  merge(.,
        dfSNPreference,
        by=c("location","nt")) %>%
  dplyr::group_by(date,sample_alias) %>%
  dplyr::mutate(n.loci.shared.per.date=length(unique(location))) %>%
  dplyr::ungroup() %>%
  dplyr::select(date,n.loci.shared.per.date,sample_alias,sample) %>%
  unique()
write.csv(tmp,"SNP_overlap.dates.csv")
