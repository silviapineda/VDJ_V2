rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis T cells from Adaptive
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: August, 2018
############################################################################################
library(gplots)
require(graphics)
library(RColorBrewer)
library(seqinr)
require(dplyr)
library(tigger)

setwd("/Users/Pinedasans/VDJ_V2/")

##Read all the files
files <- list.files("/Users/Pinedasans/VDJ_V2/Data/Adaptive_TSV/TCRB_data/")
data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("Data/Adaptive_TSV/TCRB_data/", i, sep = ""))
  t$sample<-substr(i, 1, nchar(i)-9)
  data <- rbind(data, t)
}
save(data, file="Data/TCR_adaptive.RData")


########################
#### Prepare the data ##
########################
load("Data/TCR_adaptive.RData")

##############################
### Quality Control
#############################

#dim(data)=20,413,175
data_qc<-data[which(data$cloneResolved=="VDJ"),] #20260842
##Discar STOP codons and out of frame
data_qc<-data_qc[which(data_qc$sequenceStatus=="In"),] #16380787

##Read the clinical annotations
clin_annot <- read.csv("/Users/Pinedasans/VDJ_V2/Data/ClinicalData.csv")
rownames(clin_annot) <- clin_annot$Adaptive.Sample.ID

#MERGE DATA
clin_annot$Phenotype_CADI<-factor(clin_annot$Phenotype_CADI)
data_qc$clin = clin_annot[data_qc$sample,5] ###Add the type of clinical endpoint
data_qc$clin2 = clin_annot[data_qc$sample,6] ###Add the type of clinical endpoint
data_qc$time_days = clin_annot[data_qc$sample,13] ###Add the time it was taking

data_qc$follow_up = clin_annot[data_qc$sample,3] ###Add the time it was taking

###Extract nucleotide CDR3 region
data_qc$cdr3_n<-substr(data_qc$nucleotide,data_qc$vIndex+1,data_qc$vIndex+data_qc$cdr3Length)

##Obtain the ID for the clone call
data_qc$V_J_lenghCDR3 = paste(data_qc$vGeneName,data_qc$jGeneName,data_qc$cdr3Length,sep="_")
data_qc$seq_ID<-seq(1,nrow(data_qc))
###save the data to call the clones by all samples
data_clonesInference<-data_qc[,c("seq_ID","sample","cdr3_n","cdr3Length","vGeneName","jGeneName","V_J_lenghCDR3")]
write.table(data_clonesInference,file="Data/data_for_cloneInfered_allvgenes_TCR.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides<-read.csv("Data/ClonesInfered_allvgenes_TCR.csv")
data_merge<-merge(data_qc,nucleotides[,c("seq_ID","CloneId")],by=c("seq_ID"))
##Count number of reads and  clones per sample 
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")
clones_count<- unique(data_merge[,c("sample","V_J_lenghCDR3_CloneId")])
clones<-data.matrix(table(clones_count$sample))
##Read counts per sample and data point
read_count <- table(data_merge$sample)
id_sample<-match(clin_annot$Adaptive.Sample.ID,rownames(clones))
reads_clones_annot <- cbind(clin_annot, clones[id_sample,1],read_count[id_sample])
colnames(reads_clones_annot)[c(48:50)]<-c("clones","sample","reads")
save(data_merge,reads_clones_annot,file="Data/VDJ_DataTCR.Rdata")
