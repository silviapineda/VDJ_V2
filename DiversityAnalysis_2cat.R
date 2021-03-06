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

annot_qc$Phenotype<-factor(annot_qc$Phenotype2,levels=c("NR","AR"))

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
ggplot(data=annot_qc, aes(x=time_days, y=clones, group=Sample.ID, shape=Phenotype, color=Phenotype)) +
  geom_line() +
  geom_point() +
  #ylim(0,120000) +
  scale_colour_manual(values=c("chartreuse4", "darkorange2")) +
  ggtitle("Number Clones Longitudinal")
dev.off()

# Plot for follow-up
tiff("Results/Follow_up.tiff",res=300,w=2500,h=2000)
ggplot(data=annot_qc, aes(x=followup, y=Sample.ID, shape=Phenotype, color=Phenotype)) +
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

####Association with phenotype
tiff("Results/boxplot_entropy_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0, aes(factor(Phenotype),entropy,fill=Phenotype))  + #ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of entropy")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6, aes(factor(Phenotype),entropy,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of entropy")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'entropy'
example_table<-annot_qc_time0 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_entropy_<=0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_entropy_6.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()

#############################
###Indivuduals under 18 ####
###########################
annot_qc$AgeRec<-ifelse(annot_qc$Recipient.Age..TX<=18,"<=18",">18")
annot_qc_time0_18<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$AgeRec=="<=18"),]
annot_qc_time6_18<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$AgeRec=="<=18"),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_18age.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_18, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_18, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

#############################
###Indivuduals over 18 ####
###########################
annot_qc$AgeRec<-ifelse(annot_qc$Recipient.Age..TX<=18,"<=18",">18")
annot_qc_time0_18<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$AgeRec==">18"),]
annot_qc_time6_18<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$AgeRec==">18"),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_over18age.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_18, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_18, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc_time0_18 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_0_over18.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_18 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6_over18.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()


#############################
###Indivuduals under 5 ####
###########################
annot_qc$AgeRec<-ifelse(annot_qc$Recipient.Age..TX<=5,"<=5",">5")
annot_qc_time0_5<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$AgeRec=="<=5"),]
annot_qc_time6_5<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$AgeRec=="<=5"),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6_under5age.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0_5, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6_5, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'clones'
example_table<-annot_qc_time0_18 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_0_over18.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_18 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6_over18.tiff",res=300,w=4000,h=500)
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
annot_qc_time0_24<-annot_qc[which(annot_qc$time_days<=0 & annot_qc$followup>=24),]
annot_qc_time6_24<-annot_qc[which(annot_qc$time_days>=180 & annot_qc$followup>=24),]

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
tiff("Results/Table_clones_<=0_24followip.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_24 %>% 
  finalfit(dependent, explanatory)
tiff("Results/Table_clones_6_24followup.tiff",res=300,w=4000,h=500)
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
example_table<-annot_qc_time0_3 %>% 
  finalfit(dependent, explanatory)
tiff("Results/BCR/Table_RecAge3_0.tiff",res=300,w=4000,h=500)
grid.table(example_table)
dev.off()
example_table<-annot_qc_time6_3 %>% 
  finalfit(dependent, explanatory)
tiff("Results/BCR/Table_RecAge3_6.tiff",res=300,w=4000,h=500)
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

##################################
### Match by age at transplant ###
##################################
library(e1071)
m <- matchControls(annot_qc_time0$Phenotype ~  annot_qc_time0$Recipient.Age..TX,
                                    contlabel = "NR",caselabel = "AR",replace = T)

annot_qc_matched_0<-rbind(annot_qc_time0[as.numeric(m$cases),],annot_qc_time0[as.numeric(m$controls),])

p1<-ggplot(annot_qc_matched_0, aes(factor(Phenotype),clones,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4", "darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
print(p1)

explanatory = c("Phenotype","Recipient.Age..TX")
dependent = 'clones'
example_table<-annot_qc_matched_0 %>% 
  finalfit(dependent, explanatory)
grid.table(example_table)


#################################
### Propensity score ###########
###############################
install.packages("Matching", dependencies=TRUE)
install.packages("rgenoud")

annot_qc_0<-annot_qc[which(annot_qc$time_cat=="=0"),]

library(Matching)
glm1  <- glm(Phenotype~RecipientAgeTX, family=binomial, data=annot_qc_0)
X  <- glm1$fitted
Y  <- annot_qc_0$clones
Tr  <- ifelse(annot_qc_0$Phenotype=="NR",0,1)
rr  <- Match(Y=Y, Tr=Tr, X=X, M=1);
summary(rr)
boxplot(rr$mdata$Y~rr$mdata$Tr)

####################################
##### Adaptive Analysis ###########
###################################
adaptive_data<-read.csv("Data/SampleOverview_02-12-2019_2-36-20_AM.csv")
adaptive_data$sample_name<-as.character(adaptive_data$sample_name)
adaptive_data$Adaptive.Sample.ID<-substr(adaptive_data$sample_name,1,nchar(adaptive_data$sample_name)-4)

annot_qc_adaptive<-merge(annot_qc,adaptive_data,by="Adaptive.Sample.ID")


################################
###Define time<=0 AND time>=6 ##
################################
annot_qc_time0<-annot_qc_adaptive[which(annot_qc_adaptive$time_days<=0),]
annot_qc_time6<-annot_qc_adaptive[which(annot_qc_adaptive$time_days>=180),]

####Association with phenotype
tiff("Results/boxplot_clones_0_6.tiff",h=1800,w=2000,res=300)
p4<-ggplot(annot_qc_time0, aes(factor(Phenotype),productive_rearrangements,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time <=0",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)
p5<-ggplot(annot_qc_time6, aes(factor(Phenotype),productive_rearrangements,fill=Phenotype))  + ylim(-1000,105000) +
  geom_boxplot() + scale_fill_manual(values=c("chartreuse4","darkorange2")) + #theme(text = element_text(size=15)) +
  labs(title="time >=6",x="Clinical outcome", y = "Number of clones")   + theme(legend.position="none") #+ stat_summary(fun.data=data_summary)

grid.arrange(p4,p5,ncol=3)
dev.off()

##regression
explanatory = c("Phenotype")
dependent = 'productive_rearrangements'
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

