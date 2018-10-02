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
library(lmerTest)

setwd("/Users/Pinedasans/VDJ_V2/")
load("Data/VDJ_allVgenes_DataIGH.Rdata")

###########################
## 1. Diversity measures###
###########################
# sample_unique<-unique(data_merge$sample)
# entropy<-rep(NA,length(sample_unique))
# simpson<-rep(NA,length(sample_unique))
# for (i in 1:length(sample_unique)){
#   print(i)
#   data_sample_unique<-data_merge[which(data_merge$sample==sample_unique[i]),]
#   clones_sample<-data_sample_unique[,"V_J_lenghCDR3_CloneId"]
#   #To write file to run with Recon
#   #write.delim(data.frame(table(table(clones_sample))),file=paste("Data/RECON/clones_",sample_unique[i],".txt",sep=""),sep="\t",col.names=F)
#   fi<-as.numeric(table(clones_sample))/length(clones_sample)
#   hi<-fi*log2(fi)
#   entropy[i]=-sum(hi)
#   simpson[i]=sum(fi*fi)
# }
# entropy_norm<-entropy/max(entropy,na.rm = T)
# clonality<-(1-entropy_norm)
# names(clonality)<-sample_unique
# diversity<-cbind(clonality,entropy,simpson)
# write.csv(diversity,"Data/diversity_Vgenes.csv")

diversity<-read.csv("Data/diversity_Vgenes.csv")

id<-match(reads_clones_annot$sample,diversity$X)
reads_clones_annot<-cbind(reads_clones_annot,diversity[id,2:4])

####QC on the number of clones
annot_qc<-reads_clones_annot[which(reads_clones_annot$clones>100),]
ggplot(data=annot_qc, aes(x=time_days, y=clones, group=Sample.ID, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("darkorange2", "chartreuse4")) +
  ggtitle("Number Clones Longitudinal")

###Define new variables
annot_qc$AgeRec<-ifelse(annot_qc$Recipient.Age..TX<=20,"<=20",">20")
annot_qc$AgeRec<-factor(annot_qc$AgeRec)

##barplot with ages
annot_qc_time0<-annot_qc[which(annot_qc$time_month==0),]
xx<-barplot(annot_qc_time0$Recipient.Age..TX,col = COLOR[factor(annot_qc_time0$Phenotype)],ylab = "Age Recipient")
text(x=xx,y = annot_qc_time0$Recipient.Age..TX, label = annot_qc_time0$Recipient.Age..TX, pos = 3, cex = 0.8, col = "red")
annot_qc_time6<-annot_qc[which(annot_qc$time_month>0),]
xx<-barplot(annot_qc_time6$Recipient.Age..TX,col = COLOR[factor(annot_qc_time6$Phenotype)],ylab = "Age Recipient")
text(x=xx,y = annot_qc_time6$Recipient.Age..TX, label = annot_qc_time6$Recipient.Age..TX, pos = 3, cex = 0.8, col = "red")

##################################
#####Analysis by time and clin ##
##################################
###############
###Barplots
##############
#1. clones
COLOR=c("darkorange2", "chartreuse4")
tiff("Results/barplot_clones.tiff",res=300,w=2500,h=2000)
par(mfrow=c(3,1))
barplot(annot_qc$clones[which(annot_qc$time_days<0)],
        col = COLOR[annot_qc$Phenotype[which(annot_qc$time_days<0)]],
        names.arg = annot_qc$sample[which(annot_qc$time_days<0)],
        cex.names=0.9,las=2,ylim = c(0,105000),ylab = c("Clones"))
legend(0, 105000, legend=levels(annot_qc$Phenotype[which(annot_qc$time_days>0)]),col=COLOR,pch=15, cex=1)
barplot(annot_qc$clones[which(annot_qc$time_days==0)],
        col = COLOR[annot_qc$Phenotype[which(annot_qc$time_days==0)]],
        names.arg = annot_qc$sample[which(annot_qc$time_days==0)],
        cex.names=0.9,las=2,ylim = c(0,105000),ylab = c("Clones"))
barplot(annot_qc$clones[which(annot_qc$time_days>0)],
        col = COLOR[annot_qc$Phenotype[which(annot_qc$time_days>0)]],
        names.arg = annot_qc$sample[which(annot_qc$time_days>0)],
        cex.names=0.9,las=2,ylim = c(0,105000),ylab = c("Clones"))
dev.off()


####Violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
###Boxplot
tiff("Results/boxplot_clones.tiff",h=1800,w=2000,res=300)
p1 = ggplot(annot_qc[which(annot_qc$time_days<0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days<0)]), 
                annot_qc$clones[which(annot_qc$time_days<0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + #theme(text = element_text(size=15)) +
  labs(title="time <0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p2 = ggplot(annot_qc[which(annot_qc$time_days==0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days==0)]), 
                annot_qc$clones[which(annot_qc$time_days==0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + #theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p3 = ggplot(annot_qc[which(annot_qc$time_days>0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days>0)]), 
                annot_qc$clones[which(annot_qc$time_days>0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  #theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Number of clones")  + theme(legend.position="none") #+ stat_summary(fun.data=data_summary) 

grid.arrange(p1,p2,p3,ncol=3)
dev.off()

summary(glm(annot_qc$clones[which(annot_qc$time_days<0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days<0)]))
summary(glm(annot_qc$clones[which(annot_qc$time_days==0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days==0)]))
summary(glm(annot_qc$clones[which(annot_qc$time_days>0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days>0)]))
###time<=0
tiff("Results/boxplot_clones_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc[which(annot_qc$time_month<=0),], 
       aes(factor(annot_qc$Phenotype[which(annot_qc$time_month<=0)]), 
           annot_qc$clones[which(annot_qc$time_month<=0)],fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc[which(annot_qc$time_month>=6),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_month>=6)]), 
                annot_qc$clones[which(annot_qc$time_month>=6)],fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
summary(glm(annot_qc$clones[which(annot_qc$time_month<=0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_month<=0)]))
summary(glm(annot_qc$clones[which(annot_qc$time_month>=6)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_month>=6)]))
dev.off()

#####Recipient Age
summary(glm(annot_qc$clones[which(annot_qc$time_days<=0)] 
            ~ annot_qc$Phenotype[which(annot_qc$time_days<=0)]+
              annot_qc$Recipient.Age..TX[which(annot_qc$time_days<=0)]))
summary(glm(annot_qc$clones[which(annot_qc$time_month>=6)] 
            ~ annot_qc$Phenotype[which(annot_qc$time_month>=6)]+
              annot_qc$Recipient.Age..TX[which(annot_qc$time_month>=6)]))

annot_qc_time0<-annot_qc[which(annot_qc$time_month<=0),]
annot_qc_time6<-annot_qc[which(annot_qc$time_month>=6),]

tiff("Results/plot_age.tiff",h=2000,w=2500,res=300)
par(mfrow=c(2,1))
plot(annot_qc_time0$Recipient.Age..TX,annot_qc_time0$clones,col=COLOR[annot_qc_time0$Phenotype],pch=19,
     ylab = c("Clones"),xlab = c("Age recipient"),ylim = c(0,105000),main = "time<=0")
legend(60, 104000, legend=levels(annot_qc$Phenotype),col=COLOR,pch=15, cex=1)
plot(annot_qc_time6$Recipient.Age..TX,annot_qc_time6$clones,col=COLOR[annot_qc_time6$Phenotype],pch=19,
     ylab = c("Clones"),xlab = c("Age recipient"),ylim = c(0,105000),main = "time>=6")
legend(60, 104000, legend=levels(annot_qc$Phenotype),col=COLOR,pch=15, cex=1)
dev.off()

summary(glm(annot_qc_time0$clones~annot_qc_time0$Recipient.Age..TX))
summary(glm(annot_qc_time6$clones~annot_qc_time6$Recipient.Age..TX))


###Young < 20 
annot_qc_young<-annot_qc[which(annot_qc$AgeRec=="<=20"),]
tiff("Results/boxplot_clones_young.tiff",h=1800,w=2000,res=300)
p1 = ggplot(annot_qc_young[which(annot_qc_young$time_days<0),], 
            aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days<0)]), 
                annot_qc_young$clones[which(annot_qc_young$time_days<0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + #theme(text = element_text(size=15)) +
  labs(title="time <0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p2 = ggplot(annot_qc_young[which(annot_qc_young$time_days==0),], 
            aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days==0)]), 
                annot_qc_young$clones[which(annot_qc_young$time_days==0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + #theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p3 = ggplot(annot_qc_young[which(annot_qc_young$time_days>0),], 
            aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days>0)]), 
                annot_qc_young$clones[which(annot_qc_young$time_days>0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  #theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Number of clones")  + theme(legend.position="none") #+ stat_summary(fun.data=data_summary) 

grid.arrange(p1,p2,p3,ncol=3)
dev.off()

summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days<0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days<0)]
            +annot_qc_young$Recipient.Age..TX[which(annot_qc_young$time_days<0)]))
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days==0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days==0)]
            +annot_qc_young$Recipient.Age..TX[which(annot_qc_young$time_days==0)]))
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days>0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days>0)]
            +annot_qc_young$Recipient.Age..TX[which(annot_qc_young$time_days>0)]))


#2. ENTTROPY
COLOR=c("darkorange2", "chartreuse4")
tiff("Results/barplot_entropy.tiff",res=300,w=2500,h=2000)
par(mfrow=c(2,1))
barplot(annot_qc$entropy[which(annot_qc$SampleTimeDifferencesFromTX_round.months.==0)],
        col = COLOR[annot_qc$Phenotype[which(annot_qc$SampleTimeDifferencesFromTX_round.months.==0)]],
        names.arg = annot_qc$sample[which(annot_qc$SampleTimeDifferencesFromTX_round.months.==0)],
        cex.names=0.9,las=2,ylim = c(0,14),ylab = c("entropy"))
legend(0, 14, legend=levels(annot_qc$Phenotype[which(annot_qc$SampleTimeDifferencesFromTX_round.months.>0)]),col=COLOR,pch=15, cex=1)
barplot(annot_qc$entropy[which(annot_qc$SampleTimeDifferencesFromTX_round.months.>0)],
        col = COLOR[annot_qc$Phenotype[which(annot_qc$SampleTimeDifferencesFromTX_round.months.>0)]],
        names.arg = annot_qc$sample[which(annot_qc$SampleTimeDifferencesFromTX_round.months.>0)],
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
p1 = ggplot(annot_qc[which(annot_qc$time==0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time==0)]), 
                annot_qc$entropy[which(annot_qc$time==0)],fill=Phenotype)) + ylim(0,14) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Entropy") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")

p2 = ggplot(annot_qc[which(annot_qc$time>0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time>0)]), 
                annot_qc$entropy[which(annot_qc$time>0)],fill=Phenotype)) + ylim(0,14) +
  geom_violin(adjust = .5) + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Entropy") + stat_summary(fun.data=data_summary)  + theme(legend.position="none")


grid.arrange(p1, p2,ncol=2)
dev.off()

summary(glm(annot_qc$entropy[which(annot_qc$time==0)] 
            ~ annot_qc$Phenotype[which(annot_qc$time==0)]))
summary(glm(annot_qc$entropy[which(annot_qc$time>0)] 
            ~ annot_qc$Phenotype[which(annot_qc$time>0)]))

########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
lmm <- lmer(annot_qc_young$clones ~  Phenotype*time_days 
            + (1 | Sample.ID),data=annot_qc_young)
summary(lmm)
anova(lmm)

tiff("Results/plot_lmer_clones.tiff",h=1500,w=2000,res=300)
p <- ggplot(lmm, aes(x = time_days, 
                         y = annot_qc_young$clones, colour = Phenotype)) +
  scale_colour_manual(values=c("darkorange2", "chartreuse4")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (days)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()

