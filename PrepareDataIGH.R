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
  data <- rbind(data, t)
}
save(data, file="Data/IGH_adaptive.RData")



load("Data/IGH_adaptive.RData")
