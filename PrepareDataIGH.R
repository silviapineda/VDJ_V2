rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies from Adaptive
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

setwd("/Users/Pinedasans/VDJ_V2/")

##Read all the files and convert into fasta files
files <- list.files("/Users/Pinedasans/VDJ_V2/Data/Adaptive_TSV/IGH_data/")
data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("Data/Adaptive_TSV/IGH_data/", i, sep = ""))
  seq(1,nrow(t))
  t$sample<-substr(i, 1, nchar(i)-9)
  data <- rbind(data, t)
}
save(data, file="Data/IGH_adaptive.RData")

########################
#### Prepare the data ##
########################
load("Data/IGH_adaptive.RData")

##############################
### Quality Control
#############################
##plot histogram for the V score
hist(data$V_SCORE,xlab = "V gene score", main="Histogram: Adaptive Data")
abline(v=140,col="red")

#dim(data)=10,008,758
data_qc<-data[which(data$cloneResolved=="VDJ"),] #5,874,739
##Discar STOP codons and out of frame
data_qc<-data_qc[which(data_qc$sequenceStatus=="In"),] #4,849,865

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
data_clonesInference<-data_qc[,c("SEQUENCE_ID","sample","CDR3_IGBLAST","CDR3_length","v_gene","j_gene","V_J_lenghCDR3")]
write.table(data_clonesInference,file="Data/data_for_cloneInfered_allvgenes_IGH.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides<-read.csv("Data/ClonesInfered_allvgenes_IGH.csv")
data_merge<-merge(data_qc,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))
##Count number of reads and  clones per sample 
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")
clones_count<- unique(data_merge[,c("sample","V_J_lenghCDR3_CloneId")])
clones<-data.matrix(table(clones_count$sample))
##Read counts per sample and data point
read_count <- table(data_merge$sample)
id_sample<-match(clin_annot$Adaptive.Sample.ID,rownames(clones))
reads_clones_annot <- cbind(clin_annot, clones[id_sample,1],read_count[id_sample])
colnames(reads_clones_annot)[c(48:50)]<-c("clones","sample","reads")
save(data_merge,reads_clones_annot,file="Data/VDJ_DataIGH.Rdata")
