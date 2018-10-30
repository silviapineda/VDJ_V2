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
library(finalfit)
library(dplyr)


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
  #write.delim(data.frame(table(table(clones_sample))),file=paste("Data/RECON/clones_",sample_unique[i],".txt",sep=""),sep="\t",col.names=F)
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

###########
diversity<-read.csv("Data/diversity.csv")

#####Read Recon results
#diversity_recon<-read.table("Data/test_D_number_table.txt",header=T)

#####clones
#clones_recon<-diversity_recon$est_0.0D
#entropy_recon<-log(diversity_recon$est_1.0D)
#simpson_recon<-1/diversity_recon$est_2.0D
#entropy_norm<-entropy_recon/max(entropy_recon,na.rm = T)
#clonality_recon<-(1-entropy_norm)

##Put all diversity measures together
id<-match(reads_clones_annot$Adaptive_Sample_ID,diversity$X)
reads_clones_annot<-cbind(reads_clones_annot,diversity[id,2:4])

####QC on the number of clones
annot_qc<-reads_clones_annot[which(reads_clones_annot$clones>100),]
annot_qc$Phenotype<-factor(annot_qc$Phenotype,levels = c("NoRej","bAR","AR"))
annot_qc$followup_days<-as.numeric(as.character(annot_qc$followup_days))
annot_qc$DonorAge<-as.numeric(as.character(annot_qc$DonorAge))
tiff("Results/Time_clones.tiff",res=300,w=2500,h=2000)
ggplot(data=annot_qc, aes(x=time_days, y=clones, group=sample_id, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Number Clones Longitudinal")
dev.off()

annot_qc$followup_days<-replace(annot_qc$followup_days,annot_qc$followup_days<0,0)
annot_qc$followup_month<-annot_qc$followup_days/30
tiff("Results/Follow_up.tiff",res=300,w=2500,h=2000)
ggplot(data=annot_qc, aes(x=followup_month, y=sample, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  ggtitle("Patient follow-up")
dev.off()

###Define new variables
annot_qc$AgeRec<-ifelse(annot_qc$RecipientAgeTX<=20,"<=20",">20")
annot_qc$AgeRec<-factor(annot_qc$AgeRec)
annot_qc$time_cat<-ifelse(annot_qc$time_days<0,"<0",
                          ifelse(annot_qc$time_days==0,"=0",
                                 ifelse(annot_qc$time_days>0 & annot_qc$time_days<180,"0-6",
                                        ifelse(annot_qc$time_days>=180,">=6",NA))))
annot_qc$time_cat<-factor(annot_qc$time_cat)

####Demographics####
annot_qc_unique<-annot_qc[which(duplicated(annot_qc$sample_id)==F),]
table(annot_qc_unique$Phenotype)
###follow-up time
summary(annot_qc_unique$followup_month[which(annot_qc_unique$Phenotype=="NoRej")],na.rm=T)
summary(annot_qc_unique$followup_month[which(annot_qc_unique$Phenotype=="AR")],na.rm=T)
summary(annot_qc_unique$followup_month[which(annot_qc_unique$Phenotype=="bAR")],na.rm=T)
fit = lm(annot_qc_unique$followup_month ~ annot_qc_unique$Phenotype)
anova(fit)
###Age recipient
summary(annot_qc_unique$RecipientAgeTX[which(annot_qc_unique$Phenotype=="NoRej")],na.rm=T)
summary(annot_qc_unique$RecipientAgeTX[which(annot_qc_unique$Phenotype=="AR")],na.rm=T)
summary(annot_qc_unique$RecipientAgeTX[which(annot_qc_unique$Phenotype=="bAR")],na.rm=T)
fit = lm(annot_qc_unique$RecipientAgeTX ~ annot_qc_unique$Phenotype)
anova(fit)
###Immunosupression
summary(annot_qc_unique$Immusosupression[which(annot_qc_unique$Phenotype=="NoRej")])
summary(annot_qc_unique$Immusosupression[which(annot_qc_unique$Phenotype=="AR")])
summary(annot_qc_unique$Immusosupression[which(annot_qc_unique$Phenotype=="bAR")])
chisq.test(annot_qc_unique$Immusosupression,annot_qc_unique$Phenotype)
###ESRD
summary(annot_qc_unique$ESRD[which(annot_qc_unique$Phenotype=="NoRej")],na.rm=T)
summary(annot_qc_unique$ESRD[which(annot_qc_unique$Phenotype=="AR")],na.rm=T)
summary(annot_qc_unique$ESRD[which(annot_qc_unique$Phenotype=="bAR")],na.rm=T)
chisq.test(annot_qc_unique$ESRD,annot_qc_unique$Phenotype)
###Race
summary(annot_qc_unique$RecipientRace[which(annot_qc_unique$Phenotype=="NoRej")])
summary(annot_qc_unique$RecipientRace[which(annot_qc_unique$Phenotype=="AR")])
summary(annot_qc_unique$RecipientRace[which(annot_qc_unique$Phenotype=="bAR")])
chisq.test(annot_qc_unique$RecipientRace,annot_qc_unique$Phenotype)

##DonorType
summary(annot_qc_unique$DonorType[which(annot_qc_unique$Phenotype=="NoRej")])
summary(annot_qc_unique$DonorType[which(annot_qc_unique$Phenotype=="AR")])
summary(annot_qc_unique$DonorType[which(annot_qc_unique$Phenotype=="bAR")])
chisq.test(annot_qc_unique$DonorType,annot_qc_unique$Phenotype)

##RecipientSex
summary(annot_qc_unique$RecipientSex[which(annot_qc_unique$Phenotype=="NoRej")],na.rm=T)
summary(annot_qc_unique$RecipientSex[which(annot_qc_unique$Phenotype=="AR")],na.rm=T)
summary(annot_qc_unique$RecipientSex[which(annot_qc_unique$Phenotype=="bAR")],na.rm=T)
chisq.test(annot_qc_unique$RecipientSex,annot_qc_unique$Phenotype)
chisq.test(annot_qc_unique$RecipientSex[which(annot_qc_unique$Phenotype!="bAR")],
           annot_qc_unique$Phenotype[which(annot_qc_unique$Phenotype!="bAR")])

##Time
summary(annot_qc$time_cat[which(annot_qc$Phenotype=="NoRej")],na.rm=T)
summary(annot_qc$time_cat[which(annot_qc$Phenotype=="AR")],na.rm=T)
summary(annot_qc$time_cat[which(annot_qc$Phenotype=="bAR")],na.rm=T)
chisq.test(annot_qc$time_cat,annot_qc$Phenotype)
chisq.test(annot_qc$time_cat[which(annot_qc$Phenotype!="bAR")],
           annot_qc$Phenotype[which(annot_qc$Phenotype!="bAR")])


##barplot with ages
COLOR=c("chartreuse4", "dodgerblue3","darkorange2")
annot_qc_time0<-annot_qc[which(annot_qc$time_days==0),]
tiff("Results/barplot_age_time0_Adaptive.tiff",res=300,w=3500,h=2000)
xx<-barplot(annot_qc_time0$RecipientAgeTX,col = COLOR[factor(annot_qc_time0$Phenotype)],ylab = "Age Recipient")
text(x=xx,y = annot_qc_time0$RecipientAgeTX, label = annot_qc_time0$RecipientAgeTX, pos = 3, cex = 0.8, col = "red")
dev.off()
#fit = lm(annot_qc_time0$RecipientAgeTX ~ annot_qc_time0$Phenotype)
#summary(fit)
explanatory = c("Phenotype")
dependent = 'RecipientAgeTX'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_age_time0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

tiff("Results/barplot_age_time>0_Adaptive.tiff",res=300,w=3500,h=2000)
annot_qc_time6<-annot_qc[which(annot_qc$time_days>0),]
xx<-barplot(annot_qc_time6$RecipientAgeTX,col = COLOR[factor(annot_qc_time6$Phenotype)],ylab = "Age Recipient")
text(x=xx,y = annot_qc_time6$RecipientAgeTX, label = annot_qc_time6$RecipientAgeTX, pos = 3, cex = 0.8, col = "red")
dev.off()
explanatory = c("Phenotype")
dependent = 'RecipientAgeTX'
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_age_time6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

##################################
#####Analysis by time and clin ##
##################################
###############
###Barplots
##############
#1. clones
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
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p2 = ggplot(annot_qc[which(annot_qc$time_days==0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days==0)]), 
                annot_qc$clones[which(annot_qc$time_days==0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p3 = ggplot(annot_qc[which(annot_qc$time_days>0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days>0)]), 
                annot_qc$clones[which(annot_qc$time_days>0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +  #theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Number of clones")  + theme(legend.position="none") #+ stat_summary(fun.data=data_summary) 

grid.arrange(p1,p2,p3,ncol=3)
dev.off()


summary(glm(annot_qc$clones[which(annot_qc$time_days<0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days<0)]))
summary(glm(annot_qc$clones[which(annot_qc$time_days==0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days==0)]))
summary(glm(annot_qc$clones[which(annot_qc$time_days>0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days>0)]))

explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc[which(annot_qc$time_days<=0),] %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_by_outcome.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

################################
###Define time<=0 AND time>=6 ##
################################
annot_qc_time0<-annot_qc[which(annot_qc$time_days<=0),]
annot_qc_time6<-annot_qc[which(annot_qc$time_days>=180),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

summary(glm(clones~Phenotype, data=annot_qc_time0))
summary(glm(clones~Phenotype, data=annot_qc_time6))

#############################
###Indivuduals under 20 ####
###########################
annot_qc_time0_20<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$AgeRec=="<=20"),]
annot_qc_time6_20<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$AgeRec=="<=20"),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_20age.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_20, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_20, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

summary(glm(clones~Phenotype, data=annot_qc_time0_20))
summary(glm(clones~Phenotype, data=annot_qc_time6_20))

############################################
### Association with all the variables ####
###########################################

#####Recipient Age
tiff("Results/boxplot_clones_RecAge.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(RecipientAgeTX,clones,color=Phenotype)) + 
geom_point() + geom_smooth(method='lm') + labs(title="time <=0",x="Rec.Age", y = "Clones") +
    scale_color_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
  
p2 = ggplot(annot_qc_time6,aes(RecipientAgeTX,clones,color=Phenotype)) + 
  geom_point() + geom_smooth(method='lm') + labs(title="time >=6",x="Rec.Age", y = "Clones") +
  scale_color_manual(values=c("chartreuse4", "dodgerblue3","darkorange2"))
grid.arrange(p1,p2,ncol=2)
dev.off()
summary(glm(clones~Phenotype+RecipientAgeTX, data=annot_qc_time0))
summary(glm(clones~Phenotype+RecipientAgeTX, data=annot_qc_time6))

explanatory = c("Phenotype","RecipientAgeTX")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecAge_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecAge_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#####Immunosupression
tiff("Results/boxplot_clones_Immuno.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(Immusosupression),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="Immunosuppression", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(Immusosupression),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="Immunosuppression", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

summary(glm(clones~Phenotype+Immusosupression, data=annot_qc_time0))
summary(glm(clones~Phenotype+Immusosupression, data=annot_qc_time6))

explanatory = c("Phenotype","Immusosupression")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_Immuno_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_Immuno_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()


#####DonorType
tiff("Results/boxplot_clones_DonorType.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(DonorType),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="DonorType", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(DonorType),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="DonorType", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

summary(glm(clones~Phenotype+DonorType, data=annot_qc_time0))
summary(glm(clones~Phenotype+DonorType, data=annot_qc_time6))

explanatory = c("Phenotype","DonorType")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_DonorType_0.tiff",res=300,w=4000,h=700)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_DonorType_6.tiff",res=300,w=4000,h=700)
grid.table(example_table)
dev.off()

#####RecipientSex
tiff("Results/boxplot_clones_RecipientSex.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(RecipientSex),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="RecipientSex", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(RecipientSex),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="RecipientSex", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

summary(glm(clones~Phenotype+RecipientSex, data=annot_qc_time0))
summary(glm(clones~Phenotype+RecipientSex, data=annot_qc_time6))

explanatory = c("Phenotype","RecipientSex")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecipientSex_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecipientSex_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#####RecipientRace
tiff("Results/boxplot_clones_RecipientRace.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(factor(RecipientRace),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time <=0",x="RecipientRace", y = "Clones")
p2 = ggplot(annot_qc_time6,aes(factor(RecipientRace),clones,fill=Phenotype)) + 
  geom_boxplot() + scale_fill_manual(values=COLOR)  + labs(title="time >=6",x="RecipientRace", y = "Clones")
grid.arrange(p1,p2,ncol=2)
dev.off()

summary(glm(clones~Phenotype+RecipientRace, data=annot_qc_time0))
summary(glm(clones~Phenotype+RecipientRace, data=annot_qc_time6))

explanatory = c("Phenotype","RecipientRace")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecipientRace_0.tiff",res=300,w=4000,h=900)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_RecipientRace_6.tiff",res=300,w=4000,h=900)
grid.table(example_table)
dev.off()

###Adjust by all
summary(glm(annot_qc_time0$clones~ annot_qc_time0$Phenotype + annot_qc_time0$RecipientAgeTX
            + annot_qc_time0$Immusosupression + annot_qc_time0$DonorType
            + annot_qc_time0$RecipientSex + annot_qc_time0$RecipientRace))
summary(glm(annot_qc_time6$clones~ annot_qc_time6$Phenotype + annot_qc_time6$RecipientAgeTX
            + annot_qc_time6$Immusosupression + annot_qc_time6$DonorType
            + annot_qc_time6$RecipientSex + annot_qc_time6$RecipientRace))


tiff("Results/plot_age.tiff",h=2000,w=2500,res=300)
par(mfrow=c(2,1))
plot(annot_qc_time0$RecipientAgeTX,annot_qc_time0$clones,col=COLOR[annot_qc_time0$Phenotype],pch=19,
     ylab = c("Clones"),xlab = c("Age recipient"),ylim = c(0,105000),main = "time<=0")
legend(60, 104000, legend=levels(annot_qc$Phenotype),col=COLOR,pch=15, cex=1)
plot(annot_qc_time6$RecipientAgeTX,annot_qc_time6$clones,col=COLOR[annot_qc_time6$Phenotype],pch=19,
     ylab = c("Clones"),xlab = c("Age recipient"),ylim = c(0,105000),main = "time>=6")
legend(60, 104000, legend=levels(annot_qc$Phenotype),col=COLOR,pch=15, cex=1)
dev.off()

summary(glm(annot_qc_time0$clones~annot_qc_time0$RecipientAgeTX))
summary(glm(annot_qc_time6$clones~annot_qc_time6$RecipientAgeTX))

summary(glm(annot_qc_time0$clones~annot_qc_time0$Immusosupression))
ggplot(annot_qc[which(annot_qc_time0$time_days<=0),], 
       aes(factor(annot_qc_time0$[which(annot_qc_time0$time_days<=0)]), 
           annot_qc_time0$clones[which(annot_qc_time0$time_days<=0)],fill=Immusosupression))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Immunosupression", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

###Young < 20 
annot_qc_young<-annot_qc[which(annot_qc$AgeRec=="<=20"),]
tiff("Results/boxplot_clones_young.tiff",h=1800,w=2000,res=300)
p1 = ggplot(annot_qc_young[which(annot_qc_young$time_days<0),], 
            aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days<0)]), 
                annot_qc_young$clones[which(annot_qc_young$time_days<0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p2 = ggplot(annot_qc_young[which(annot_qc_young$time_days==0),], 
            aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days==0)]), 
                annot_qc_young$clones[which(annot_qc_young$time_days==0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p3 = ggplot(annot_qc_young[which(annot_qc_young$time_days>0),], 
            aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days>0)]), 
                annot_qc_young$clones[which(annot_qc_young$time_days>0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +  #theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Number of clones")  + theme(legend.position="none") #+ stat_summary(fun.data=data_summary) 

grid.arrange(p1,p2,p3,ncol=3)
dev.off()

summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days<0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days<0)]))
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days==0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days==0)]))
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days>0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days>0)]))

###time<=0 AND time>=6
tiff("Results/boxplot_young_clones_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_young[which(annot_qc_young$time_days<=0),], 
           aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days<=0)]), 
               annot_qc_young$clones[which(annot_qc_young$time_days<=0)],fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_young[which(annot_qc_young$time_days>=180),], 
           aes(factor(annot_qc_young$Phenotype[which(annot_qc_young$time_days>=180)]), 
               annot_qc_young$clones[which(annot_qc_young$time_days>=180)],fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days<=0)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days<=0)]))
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days>=180)] ~ 
              annot_qc_young$Phenotype[which(annot_qc_young$time_days>=180)]))
dev.off()

#####Recipient Age
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_days<=0)] 
            ~ annot_qc_young$Phenotype[which(annot_qc_young$time_days<=0)]+
              annot_qc_young$Recipient.Age..TX[which(annot_qc_young$time_days<=0)]))
summary(glm(annot_qc_young$clones[which(annot_qc_young$time_month>=6)] 
            ~ annot_qc_young$Phenotype[which(annot_qc_young$time_month>=6)]+
              annot_qc_young$Recipient.Age..TX[which(annot_qc_young$time_month>=6)]))

#2. ENTTROPY
COLOR=c("darkorange2", "chartreuse4")
###Boxplot
tiff("Results/boxplot_entropy.tiff",h=1800,w=2000,res=300)
p1 = ggplot(annot_qc[which(annot_qc$time_days<0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days<0)]), 
                annot_qc$entropy[which(annot_qc$time_days<0)],fill=Phenotype)) + ylim(5,17) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + 
  labs(title="time <0",x="Clinical outcome", y = "Number of entropy")   + theme(legend.position="none") 
p2 = ggplot(annot_qc[which(annot_qc$time_days==0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days==0)]), 
                annot_qc$entropy[which(annot_qc$time_days==0)],fill=Phenotype)) + ylim(5,17) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + 
  labs(title="time 0",x="Clinical outcome", y = "Number of entropy")   + theme(legend.position="none")
p3 = ggplot(annot_qc[which(annot_qc$time_days>0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days>0)]), 
                annot_qc$entropy[which(annot_qc$time_days>0)],fill=Phenotype)) + ylim(5,17) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) +  
  labs(title="time >0",x="Clinical outcome", y = "Number of entropy")  + theme(legend.position="none")

grid.arrange(p1,p2,p3,ncol=3)
dev.off()

summary(glm(annot_qc$entropy[which(annot_qc$time_days<0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days<0)]))
summary(glm(annot_qc$entropy[which(annot_qc$time_days==0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days==0)]))
summary(glm(annot_qc$entropy[which(annot_qc$time_days>0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_days>0)]))
###time<=0 AND time>=6
tiff("Results/boxplot_entropy_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc[which(annot_qc$time_month<=0),], 
           aes(factor(annot_qc$Phenotype[which(annot_qc$time_month<=0)]), 
               annot_qc$entropy[which(annot_qc$time_month<=0)],fill=Phenotype))  + ylim(5,17) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + 
  labs(title="time <=0",x="Clinical outcome", y = "Number of entropy")   + theme(legend.position="none") 
p5<-ggplot(annot_qc[which(annot_qc$time_month>=6),], 
           aes(factor(annot_qc$Phenotype[which(annot_qc$time_month>=6)]), 
               annot_qc$entropy[which(annot_qc$time_month>=6)],fill=Phenotype))  + ylim(5,17) +
  geom_boxplot() + scale_fill_manual(values=c("darkorange2", "chartreuse4")) + 
  labs(title="time >=6",x="Clinical outcome", y = "Number of entropy")   + theme(legend.position="none") 

grid.arrange(p4,p5,ncol=3)
summary(glm(annot_qc$entropy[which(annot_qc$time_month<=0)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_month<=0)]))
summary(glm(annot_qc$entropy[which(annot_qc$time_month>=6)] ~ 
              annot_qc$Phenotype[which(annot_qc$time_month>=6)]))
dev.off()

#####Recipient Age
summary(glm(annot_qc$entropy[which(annot_qc$time_days<=0)] 
            ~ annot_qc$Phenotype[which(annot_qc$time_days<=0)]+
              annot_qc$Recipient.Age..TX[which(annot_qc$time_days<=0)]))
summary(glm(annot_qc$entropy[which(annot_qc$time_month>=6)] 
            ~ annot_qc$Phenotype[which(annot_qc$time_month>=6)]+
              annot_qc$Recipient.Age..TX[which(annot_qc$time_month>=6)]))

########################################
###  Fitthing a longitudinal model  ####
#######################################

##Clones
lmm <- lmer(clones ~  Phenotype*time_days 
            + (1 | sample_id),data=annot_qc)
summary(lmm)
anova(lmm)

tiff("Results/plot_lmer_clones.tiff",h=1500,w=2000,res=300)
p <- ggplot(lmm, aes(x = time_days, 
                         y = annot_qc$clones, colour = Phenotype)) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue3","darkorange2")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (days)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()

##Clones young 
lmm <- lmer(clones ~  Phenotype*time_days 
            + (1 | Sample.ID),data=annot_qc_young)
summary(lmm)
anova(lmm)

tiff("Results/plot_lmer_clones_young.tiff",h=1500,w=2000,res=300)
p <- ggplot(lmm, aes(x = time_days, 
                     y = annot_qc_young$clones, colour = Phenotype)) +
  scale_colour_manual(values=c("darkorange2", "chartreuse4")) +
  geom_point(size=1.2) +
  geom_smooth(method="lm",size=0.8) +
  labs(x = "Time (days)",y = "Number of clones") + theme_bw() + theme_light()

print(p)
dev.off()
