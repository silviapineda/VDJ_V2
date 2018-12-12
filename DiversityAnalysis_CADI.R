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

COLOR=c("chartreuse4","darkorange2")


setwd("/Users/Pinedasans/VDJ_V2/")
load("Data/annot_qc.Rdata")
annot_qc$Phenotype<-factor(annot_qc$Phenotype2,levels = c("NR","AR"))

# Table 1 - Patient demographics by variable of interest ----
annot_qc_unique<-annot_qc[which(duplicated(annot_qc$Sample.ID)==F),]
table(annot_qc_unique$Phenotype)
explanatory = c("RecipientAgeTX", "Immusosupression..1.SF..2.SB.","Donor.Type", "Recipient.Sex","followup")
dependent = 'Phenotype'
table<-annot_qc_unique %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
tiff("Results/Table_demographic.tiff",res=300,w=4000,h=1200)
grid.table(table)
dev.off()

# Table 2 - Time by variable of interest ----
explanatory = c("time_cat")
dependent = 'Phenotype'
table<-annot_qc %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
tiff("Results/Table_time.tiff",res=300,w=4000,h=1200)
grid.table(table)
dev.off()

# Table 3 - Patient demographics by variable of interest and time point <=0
annot_qc_unique_0<-annot_qc[which(annot_qc$time_days<=0),]
explanatory = c("RecipientAgeTX", "Immusosupression..1.SF..2.SB.","Donor.Type", "Recipient.Sex","followup")
dependent = 'Phenotype'
table<-annot_qc_unique_0 %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
tiff("Results/Table_demographic_<=0.tiff",res=300,w=4000,h=1200)
grid.table(table)
dev.off()

# Table 4 - Patient demographics by variable of interest and time point >=6
annot_qc_unique_6<-annot_qc[which(annot_qc$time_days>=180),]
explanatory = c("RecipientAgeTX", "Immusosupression..1.SF..2.SB.","Donor.Type", "Recipient.Sex","followup")
dependent = 'Phenotype'
table<-annot_qc_unique_6 %>%
  summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)
tiff("Results/Table_demographic_>=6.tiff",res=300,w=4000,h=1200)
grid.table(table)
dev.off()

# Plot for time
tiff("Results/Time_clones.tiff",res=300,w=2500,h=2000)
ggplot(data=annot_qc, aes(x=time_days, y=clones, group=sample_id, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "darkorange2")) +
  ggtitle("Number Clones Longitudinal")
dev.off()

# Plot for follow-up
annot_qc$followup_days<-replace(annot_qc$followup_days,annot_qc$followup_days<0,0)
annot_qc$followup_month<-annot_qc$followup_days/30
tiff("Results/Follow_up.tiff",res=300,w=2500,h=2000)
ggplot(data=annot_qc, aes(x=followup_month, y=sample, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "darkorange2")) +
  ggtitle("Patient follow-up")
dev.off()

##################################
#####Analysis by time and clin ##
##################################
##Boxplot
tiff("Results/boxplot_clones.tiff",h=1800,w=2000,res=300)
p1 = ggplot(annot_qc[which(annot_qc$time_days<0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days<0)]), 
                annot_qc$clones[which(annot_qc$time_days<0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p2 = ggplot(annot_qc[which(annot_qc$time_days==0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days==0)]), 
                annot_qc$clones[which(annot_qc$time_days==0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time 0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p3 = ggplot(annot_qc[which(annot_qc$time_days>0),], 
            aes(factor(annot_qc$Phenotype[which(annot_qc$time_days>0)]), 
                annot_qc$clones[which(annot_qc$time_days>0)],fill=Phenotype)) + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) +  #theme(text = element_text(size=15)) +
  labs(title="time >0",x="Clinical outcome", y = "Number of clones")  + theme(legend.position="none") #+ stat_summary(fun.data=data_summary) 

grid.arrange(p1,p2,p3,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc[which(annot_qc$time_days<0),] %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_<0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc[which(annot_qc$time_days==0),] %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc[which(annot_qc$time_days>0),] %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_>0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

################################
###Define time<=0 AND time>=6 ##
################################
annot_qc_time0<-annot_qc[which(annot_qc$time_days<=0),]
annot_qc_time6<-annot_qc[which(annot_qc$time_days>=180),]
summary(glm(annot_qc_time0$clones ~ annot_qc_time0$Phenotype))

####Association with phenotype
tiff("Results/boxplot_clones_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_<=0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#############################
###Indivuduals under 20 ####
###########################
annot_qc_time0_20<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$AgeRec=="<=20"),]
annot_qc_time6_20<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$AgeRec=="<=20"),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_20age.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_20, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_20, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc_time0_20 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_<=0_20.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_20 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6_20.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#############################
###Indivuduals over 3  ####
###########################
annot_qc_time0_3<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$RecipientAgeTX>3),]
annot_qc_time6_3<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$RecipientAgeTX>3),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_3age.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_3, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_3, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype","RecipientAgeTX")
dependent = 'clones'
example_table<-annot_qc_time0_3 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_<=0_3.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_3 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6_3.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#######################
###Long follow-up ####
######################
annot_qc_time0_24<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$followup_month>=6),]
annot_qc_time6_24<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$followup_month>=6),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_24followup.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_24, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_24, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc_time0_24 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_<=0_20.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_24 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6_20.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()


############################################
### Association with all the variables ####
###########################################

#####Recipient Age
tiff("Results/boxplot_clones_RecAge.tiff",h=2000,w=3000,res=300)
p1 = ggplot(annot_qc_time0,aes(RecipientAgeTX,clones,color=Phenotype)) + 
  geom_point() + geom_smooth(method='lm') + labs(title="time <=0",x="Rec.Age", y = "Clones") +
  scale_color_manual(values=c("chartreuse4","darkorange2"))

p2 = ggplot(annot_qc_time6,aes(RecipientAgeTX,clones,color=Phenotype)) + 
  geom_point() + geom_smooth(method='lm') + labs(title="time >=6",x="Rec.Age", y = "Clones") +
  scale_color_manual(values=c("chartreuse4", "darkorange2"))
grid.arrange(p1,p2,ncol=2)
dev.off()

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

##########Considering all variables
explanatory = c("Phenotype","RecipientAgeTX", "Immusosupression","DonorType", "RecipientSex")
dependent = 'clones'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_All_0.tiff",res=300,w=4000,h=1000)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_All_6.tiff",res=300,w=4000,h=1000)
grid.table(example_table)
dev.off()

########## Meta-Analaysis######
library(metafor)
yi<-c(1592.3,5393)
sei<-c(479.7,3194)
res <- rma(yi,sei=sei,method="DL")
forest(res)



