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
files <- list.files("/Users/Pinedasans/VDJ_V2/Data/Adaptive/IGH_data/")
data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("Data/Adaptive/IGH_data/", i, sep = ""))
  seq(1,nrow(t))
  t$Sample <- paste0(substr(i, 1, nchar(i)-8),"_",seq(1,nrow(t)))
  t_vdj<-t[which(t$cloneResolved=="VDJ"),]
  t_vdjData = dplyr::data_frame(name = t_vdj$Sample,
                                  seq = as.character(t_vdj$nucleotide))
  writeFasta(t_vdjData, paste0("Data/",substr(i, 1, nchar(i)-8),".fasta"))
  #data <- rbind(data, t)
}

##Read the files generated with IgBlast and save into a Rdata
files <- list.files("/Users/Pinedasans/VDJ_V2/Data/IGBLAST/")
data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("Data/IGBLAST/", i, sep = ""))
  t$sample<-substr(i, 1, nchar(i)-18)
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

##Discard all the V_score < 140
data_qc<-data[which(data$V_SCORE>=140),]
##Discard the non-functional sequences
data_qc<-data_qc[which(data_qc$FUNCTIONAL=="TRUE"),]

##Read the clinical annotations
clin_annot <- read.csv("/Users/Pinedasans/VDJ_V2//Data/ClinicalData.csv")
rownames(clin_annot) <- clin_annot$Adaptive.Sample.ID

#MERGE DATA
clin_annot$Phenotype<-factor(clin_annot$Phenotype)
data_qc$clin = clin_annot[data_qc$sample,3] ###Add the type of clinical endpoint
data_qc$time_months = clin_annot[data_qc$sample,9] ###Add the time it was taking
data_qc$time_months_round = clin_annot[data_qc$sample,10] ###Add the time it was taking

##Extract the gene from the segment with the allele
data_qc$v_gene <- gsub("\\*", "", substr(data_qc$V_CALL, 1, 8))
data_qc$j_gene <- gsub("\\*", "", substr(data_qc$J_CALL, 1, 5))
data_qc$d_gene <- gsub("\\*", "", substr(data_qc$D_CALL, 1, 8))

##count the CDR3 length
data_qc$CDR3_length<-nchar(as.character(data_qc$CDR3_IGBLAST))
data_qc<-data_qc[which(data_qc$CDR3_length>0),]
##Obtain the ID for the clone call
data_qc$V_J_lenghCDR3 = paste(data_qc$v_gene,data_qc$j_gene,data_qc$CDR3_length,sep="_")

###save the data to call the clones by all samples
data_clonesInference<-data_qc[,c("SEQUENCE_ID","sample","CDR3_IGBLAST","CDR3_length","v_gene","j_gene","V_J_lenghCDR3")]
write.table(data_clonesInference,file="Data/data_for_cloneInfered_IGH.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides<-read.csv("Data/ClonesInfered_IGH.csv")
data_merge<-merge(data_qc,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))
##Count number of reads and  clones per sample 
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")
clones_count<- unique(data_merge[,c("sample","V_J_lenghCDR3_CloneId")])
clones<-data.matrix(table(clones_count$sample))
##Read counts per sample and data point
read_count <- table(data_merge$sample)
id_sample<-match(clin_annot$Adaptive.Sample.ID,rownames(clones))
reads_clones_annot <- cbind(clin_annot, clones[id_sample,1],read_count[id_sample])
colnames(reads_clones_annot)[44:46]<-c("clones","sample","reads")
save(data_merge,reads_clones_annot,file="Data/VDJ_DataIGH.Rdata")
