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
load("Data/VDJ_allVgenes_DataIGH.Rdata")

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
###Convert into 0 all those who has a low expression (<0.05)
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
p.value<-matrix(NA,3,ncol(vusage_filter))
for(i in 1:ncol(vusage_filter)){
  #time0
  fit0<-summary(glm(vusage_filter[which(vgenes$time_month<=0),i] ~ vgenes$phenotype[which(vgenes$time_month<=0)]))
  p.value[1,i]<-fit0$coefficients[2,4]
  #time>0
  fit0<-summary(glm(vusage_filter[which(vgenes$time_month>0),i] ~ vgenes$phenotype[which(vgenes$time_month>0)]))
  p.value[2,i]<-fit0$coefficients[2,4]
  #all
  fit0<-summary(glm(vusage_filter[,i] ~ vgenes$phenotype))
  p.value[3,i]<-fit0$coefficients[2,4]
}


colnames(p.value)<-colnames(vusage_filter)
rownames(p.value)<-c("time0","time>0","all")
p.value<0.05

tiff("Results/boxplot_IGHV3-23.tiff",res=300,w=1500,h=1000)
par(mfrow=c(1,3))
boxplot(vusage_filter[which(vgenes$time_month<=0),"IGHV3-23"]~vgenes$phenotype[which(vgenes$time_month<=0)],
        col=c("darkorange2","chartreuse4"),ylim=c(0.0,0.2),ylab="IHGV3-23 expression",main="time <= 0")
boxplot(vusage_filter[which(vgenes$time_month>0),"IGHV3-23"]~vgenes$phenotype[which(vgenes$time_month>0)],
        col=c("darkorange2","chartreuse4"),ylim=c(0.0,0.2),ylab="IHGV3-23 expression",main="time >0")
boxplot(vusage_filter[,"IGHV3-23"]~vgenes$phenotype,
        col=c("darkorange2","chartreuse4"),ylim=c(0.0,0.2),ylab="IHGV3-23 expression",main="all")
dev.off()

###################
##Plot Heatmap ###
##################
 
vgenes_sign_0<-names(which(p.value[1,]<0.06))
vgenes_sign_6<-names(which(p.value[2,]<0.06))
vgenes_sign_all<-names(which(p.value[3,]<0.06))

##Heatmap with significant results
id_sign_0<-match(vgenes_sign_0,colnames(vusage_filter))
vusage_sign_0<-vusage_filter[,id_sign_0]
id_sign_6<-match(vgenes_sign_6,colnames(vusage_filter))
vusage_sign_6<-vusage_filter[,id_sign_6]
id_sign_all<-match(vgenes_sign_all,colnames(vusage_filter))
vusage_sign_all<-vusage_filter[,id_sign_all]

vusage_sign_0<-vusage_sign_0[which(vgenes$time_month<=0),]
vusage_sign_6<-vusage_sign_6[which(vgenes$time_month>0),]

annotation_col_0 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]),
  source = factor(vgenes_filter$Source[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]))
annotation_col_6 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]),
  source = factor(vgenes_filter$Source[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]))
annotation_col_24 = data.frame(
  clin = factor(vgenes_filter$clin[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  ESRD = factor(vgenes_filter$ESRD[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  immunosuppression = factor(vgenes_filter$immunosuppression[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]),
  source = factor(vgenes_filter$Source[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]))


rownames(annotation_col_0)<-vgenes_filter$Individual_id[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]
rownames(annotation_col_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]
rownames(annotation_col_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]

rownames(vusage_sign_0)<-vgenes_filter$Individual_id[which(vgenes_filter$time==0 & vgenes_filter$clin!="PNR")]
rownames(vusage_sign_6)<-vgenes_filter$Individual_id[which(vgenes_filter$time==6 & vgenes_filter$clin!="PNR")]
rownames(vusage_sign_24)<-vgenes_filter$Individual_id[which(vgenes_filter$time==24 & vgenes_filter$clin!="PNR")]



fill=c("chartreuse4","darkorange2")
fill2=brewer.pal(3,"Set2")
#fill3 = brewer.pal(3,"Set1")
#fill4 = b=brewer.pal(3,"Set3")

ann_colors = list (clin = c("NP" = fill[1], "PR" = fill[2]),
                   ESRD = c("Other/Uk" = fill2[1], "Reflux" = fill2[2], "Non-Immune/Strutural" = fill2[3]))
#immunosuppression = c("Steroid-free" = fill3[1], "Steroid-based" = fill3[2]),
#source = c("Cadaver" = fill4[1],"Living/related"=fill4[2]))

colfunc<-colorRampPalette(c("white","red4"))
tiff("heatmap_vgene_0.tiff",res=300,w=1500,h=1200)
pheatmap(t(vusage_sign_0),annotation_col = annotation_col_0,fontsize = 8,main="time 0"
         ,annotation_colors = ann_colors,rownames = F,color =  colfunc(100),border_color=F)
dev.off()

