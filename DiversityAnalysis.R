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
library(entropy)
library(ggplot2)
library(untb)
library(lme4)
library(caroline)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(lattice)


setwd("/Users/Pinedasans/VDJ_V2/")
load("Data/VDJ_DataIGH.Rdata")
###########################
## 1. Diversity measures###
###########################
sample_unique<-unique(data_merge$sample)
entropy<-rep(NA,length(sample_unique))
simpson<-rep(NA,length(sample_unique))
for (i in 1:length(sample_unique)){
  print(i)
  data_sample_unique<-data_merge[which(data_merge$sample==sample_unique[i]),]
  clones_sample<-data_sample_unique[,"V_J_lenghCDR3_CloneId"]
  #To write file to run with Recon
  write.delim(data.frame(table(table(clones_sample))),file=paste("Data/RECON/clones_",sample_unique[i],".txt",sep=""),sep="\t",col.names=F)
  fi<-as.numeric(table(clones_sample))/length(clones_sample)
  hi<-fi*log2(fi)
  entropy[i]=-sum(hi) 
  simpson[i]=sum(fi*fi) 
}
entropy_norm<-entropy/max(entropy,na.rm = T)
clonality<-(1-entropy_norm)
names(clonality)<-sample_unique
diversity<-cbind(clonality,entropy,simpson)
write.csv(diversity,"Data/diversity.csv")

id<-match(reads_clones_annot$sample,rownames(diversity))
reads_clones_annot<-cbind(reads_clones_annot,diversity[id,])

####QC on the number of reads
reads_clones_annot_reads100<-reads_clones_annot[which(reads_clones_annot$reads>100),]
ggplot(data=reads_clones_annot_reads100, aes(x=SampleTimeDifferencesFromTX_round.months., y=clones, group=sample, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("darkorange2", "chartreuse4")) +
  ggtitle("Number Clones Longitudinal")

##################################
#####Analysis by time and clin ##
##################################
###############
###Barplots
##############
#1. clones
COLOR=c("darkorange2", "chartreuse4")
tiff("Results/barplot_clones.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(reads_clones_annot_reads100$clones[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)],
        col = COLOR[reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)]],
        names.arg = reads_clones_annot_reads100$sample[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)],
        cex.names=0.9,las=2,ylim = c(0,15000),ylab = c("Clones"))
legend(0, 20000, legend=levels(reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)]),col=COLOR,pch=15, cex=1)
barplot(reads_clones_annot_reads100$clones[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)],
        col = COLOR[reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)]],
        names.arg = reads_clones_annot_reads100$sample[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)],
        cex.names=0.9,las=2,ylim = c(0,15000),ylab = c("Clones"))
dev.off()


####Violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
tiff("Results/violin_clones.tiff",h=1800,w=3000,res=300)
p1 = ggplot(reads_clones_annot_reads100[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0),], 
            aes(factor(reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)]), 
                reads_clones_annot_reads100$clones[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)],fill=Phenotype)) + ylim(0,15000) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")

p2 = ggplot(reads_clones_annot_reads100[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0),], 
            aes(factor(reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)]), 
                reads_clones_annot_reads100$clones[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)],fill=Phenotype)) + ylim(0,15000) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Number of clones") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")


grid.arrange(p1, p2,ncol=2)
dev.off()

summary(glm(reads_clones_annot_reads100$clones[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)] 
            ~ reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)]))
summary(glm(reads_clones_annot_reads100$clones[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==6)] 
            ~ reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==6)]))

table(reads_clones_annot_reads100$Phenotype,reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.)

#2. ENTTROPY
COLOR=c("darkorange2", "chartreuse4")
tiff("Results/barplot_entropy.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(reads_clones_annot_reads100$entropy[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)],
        col = COLOR[reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)]],
        names.arg = reads_clones_annot_reads100$sample[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)],
        cex.names=0.9,las=2,ylim = c(0,14),ylab = c("entropy"))
legend(0, 14, legend=levels(reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)]),col=COLOR,pch=15, cex=1)
barplot(reads_clones_annot_reads100$entropy[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)],
        col = COLOR[reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)]],
        names.arg = reads_clones_annot_reads100$sample[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)],
        cex.names=0.9,las=2,ylim = c(0,14),ylab = c("entropy"))
dev.off()


####Violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
tiff("Results/violin_entropy.tiff",h=1800,w=3000,res=300)
p1 = ggplot(reads_clones_annot_reads100[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0),], 
            aes(factor(reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)]), 
                reads_clones_annot_reads100$entropy[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)],fill=Phenotype)) + ylim(0,14) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Entropy") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")

p2 = ggplot(reads_clones_annot_reads100[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0),], 
            aes(factor(reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)]), 
                reads_clones_annot_reads100$entropy[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.>0)],fill=Phenotype)) + ylim(0,14) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Entropy") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")


grid.arrange(p1, p2,ncol=2)
dev.off()

summary(glm(reads_clones_annot_reads100$entropy[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)] 
            ~ reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==0)]))
summary(glm(reads_clones_annot_reads100$entropy[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==6)] 
            ~ reads_clones_annot_reads100$Phenotype[which(reads_clones_annot_reads100$SampleTimeDifferencesFromTX_round.months.==6)]))

########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
fm_null <- lmer(reads_clones_annot_reads100$entropy ~ Phenotype + SampleTimeDifferencesFromTX_round.months. 
                + (1 | Sample.ID),data=reads_clones_annot_reads100,REML = F)
fm_full <- lmer(reads_clones_annot_reads100$entropy ~  Phenotype*SampleTimeDifferencesFromTX_round.months. 
                + (1 | Sample.ID) ,data=reads_clones_annot_reads100,REML = F)
anova(fm_full, fm_null) 

tiff("Results/plot_lmer_entropy.tiff",h=1200,w=1400,res=300)
p <- ggplot(fm_full, aes(x = SampleTimeDifferencesFromTX_round.months., 
                         y = reads_clones_annot_reads100$entropy, colour = Phenotype)) +
  scale_colour_manual(values=c("darkorange2", "chartreuse4")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (months)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()

