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
### Date: September, 2018
############################################################################################
library(circlize)
library("RColorBrewer")
library(gtools)
library(lme4)
library("randomForest")
library("VSURF")
library(pheatmap)
library(ggplot2)
library(reshape2)

setwd("/Users/Pinedasans/VDJ_V2/")
load("Data/VDJ_DataIGH.Rdata")

#For this analysis we are putting a cut-off on clones because v-genes can be biased at low clonality 
reads_clones_annot_clones100<-reads_clones_annot[which(reads_clones_annot$clones>100),]

id<-match(data_merge$sample,reads_clones_annot_clones100$sample)
data_merge_qc<-data_merge[which(is.na(id)==F),]
data_merge_qc$sample<-factor(data_merge_qc$sample)

vgenes<-as.data.frame(unclass(table(data_merge_qc$sample,data_merge_qc$v_gene)))
id.spec<-match(rownames(vgenes),reads_clones_annot_clones100$sample)
vgenes<-cbind(vgenes,reads_clones_annot_clones100$Phenotype[id.spec],reads_clones_annot_clones100$time_days[id.spec],
              reads_clones_annot_clones100$time_month[id.spec],reads_clones_annot_clones100$Sample.ID[id.spec],
              reads_clones_annot_clones100$sample[id.spec],reads_clones_annot_clones100$Recipient.Age..TX[id.spec])
colnames(vgenes)[71:76]<-c("phenotype","time_days","time_month","Sample_id","Sample","Age")

vusage<-matrix(NA,nrow(vgenes),(ncol(vgenes)-6))
for (i in 1:70){
  vusage[,i]<-vgenes[,i]/reads_clones_annot_clones100$clones[id.spec]
}
colnames(vusage)<-colnames(vgenes)[1:70]
rownames(vusage)<-vgenes$Sample

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
xx<-replace(vusage,vusage<0.05,0)
###Those who are in lesss than 10%
vusage_filter<-vusage[,which(apply(xx,2,function(x) sum(x==0))<=267)]
##10 genes in total

AMR_usage<-colSums(vusage_filter[which(vgenes$phenotype=="AMR"),])
NoRej_usage<-colSums(vusage_filter[which(vgenes$phenotype=="NoRej"),])
AMR_usage[order(AMR_usage)]
NoRej_usage[order(NoRej_usage)]

#########
## Statistical analysis to find significant vgenes
########

vusage_filter<-vusage_filter[which(vgenes$Age<=20),]
vgenes<-vgenes[which(vgenes$Age<=20),]
##ANOVA by time points
p.value<-matrix(NA,3,10)
for(i in 1:10){
  #time0
  fit0<-summary(glm(vusage_filter[which(vgenes$time_days==0),i] ~ vgenes$phenotype[which(vgenes$time_days==0)]))
  p.value[1,i]<-fit0$coefficients[2,4]
  #time>0
  fit0<-summary(glm(vusage_filter[which(vgenes$time_days>0),i] ~ vgenes$phenotype[which(vgenes$time_days>0)]))
  p.value[2,i]<-fit0$coefficients[2,4]
  #all
  fit0<-summary(glm(vusage_filter[,i] ~ vgenes$phenotype))
  p.value[3,i]<-fit0$coefficients[2,4]
}

colnames(p.value)<-colnames(vusage_filter)
rownames(p.value)<-c("time0","time>0","all")
p.value<0.05

tiff("Results/boxplot_IGHV3-23_young.tiff",res=300,w=1500,h=1000)
par(mfrow=c(1,3))
boxplot(vusage_filter[which(vgenes$time_days==0),"IGHV3-23"]~vgenes$phenotype[which(vgenes$time_days==0)],
        col=c("darkorange2","chartreuse4"),ylim=c(0.0,0.2),ylab="IHGV3-23 expression",main="time 0")
boxplot(vusage_filter[which(vgenes$time_days>00),"IGHV3-23"]~vgenes$phenotype[which(vgenes$time_days>0)],
        col=c("darkorange2","chartreuse4"),ylim=c(0.0,0.2),ylab="IHGV3-23 expression",main="time >0")
boxplot(vusage_filter[,"IGHV3-23"]~vgenes$phenotype,
        col=c("darkorange2","chartreuse4"),ylim=c(0.0,0.2),ylab="IHGV3-23 expression",main="all")
dev.off()

 

