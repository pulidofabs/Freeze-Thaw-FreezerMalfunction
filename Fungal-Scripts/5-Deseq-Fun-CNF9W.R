#Reset R's Brain
rm(list=ls())

#Set working directory................................................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/NoCNF9WT2/")

#Load Packages-------------------------------------------------------------------------------
library("qiime2R")
library(DESeq2)
library(ggplot2)#Plotting
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(ochRe)# color paletteibrary(multcomp)#multiple comparisons--for glm 
library(tidyverse)
library(tidyr)

Metadata<-read_tsv("Fungi/Metadata/MetaRareFun-CNF9WT2.tsv")#Rare metadata
Table<-read_qza("Qiime/Merged/Fungi/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")
Tree<-read_qza("Qiime/Merged/Fungi/CoreMetrics/Phylo/Fungi-rooted-tree-1-4RT-9WT2.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Merged/Fungi/Taxonomy/Fungi-taxonomy-1-4RT-9WT2-Ver8.99.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                                                     c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#Create phyloseq artifact.........................................................................
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SeqSampleID")))





#----
#----
##----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("k__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])


sample_data(physeq)$SampleType<-factor(sample_data(physeq)$SampleType) # make as factor

#Subset data to look at sample types independently ----------------------------------------------------------------------
physeqO<-subset_samples(physeq,SampleType=="Original");physeqO#2607x19
physeqT<-subset_samples(physeq,SampleType=="Thawed");physeqT#2607x19
physeqF<-subset_samples(physeq,SampleType=="Frozen");physeqF#2607x19


#Subset each sampletype by timepoint, start with otirinals..........................................
unique(sample_data(physeqO)$Timepoint)#"T3" "T4" "T2" "T8" "T5"

#Original
physeqO2<-subset_samples(physeqO,Timepoint=="T2");sample_names(physeqO2)#11595x31
physeqO3<-subset_samples(physeqO,Timepoint=="T3");sample_names(physeqO3)#11595x36
physeqO4<-subset_samples(physeqO,Timepoint=="T4");sample_names(physeqO4)#11595x36
physeqO5<-subset_samples(physeqO,Timepoint=="T5");sample_names(physeqO5)#11595x36
physeqO8<-subset_samples(physeqO,Timepoint=="T8");sample_names(physeqO8)#11595x36


#Thawed
physeqT2<-subset_samples(physeqT,Timepoint=="T2");sample_names(physeqT2)#11595x31
physeqT3<-subset_samples(physeqT,Timepoint=="T3");sample_names(physeqT3)#11595x36
physeqT4<-subset_samples(physeqT,Timepoint=="T4");sample_names(physeqT4)#11595x36
physeqT5<-subset_samples(physeqT,Timepoint=="T5");sample_names(physeqT5)#11595x36
physeqT8<-subset_samples(physeqT,Timepoint=="T8");sample_names(physeqT8)#11595x36

#Frozen
physeqF2<-subset_samples(physeqF,Timepoint=="T2");sample_names(physeqF2)#11595x31
physeqF3<-subset_samples(physeqF,Timepoint=="T3");sample_names(physeqF3)#11595x36
physeqF4<-subset_samples(physeqF,Timepoint=="T4");sample_names(physeqF4)#11595x36
physeqF5<-subset_samples(physeqF,Timepoint=="T5");sample_names(physeqF5)#11595x36
physeqF8<-subset_samples(physeqF,Timepoint=="T8");sample_names(physeqF8)#11595x36



#----
#-----
#****************************************************************************************************-----
#--- RUN DESEQ FOR ORIGINAL SAMPLES
#****************************************************************************************************-----

#Timepoint 1 **********************************************************************----
head(sample_data(physeqO2)$Treatment, 100)

#Ran in case there are any n/a in data
physeqO2<-subset_samples(physeqO2, Treatment!= "None")

#Import phyloseq data to create a deseq object.....................
dds1<-phyloseq_to_deseq2(physeqO2, ~ Treatment)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans1 = apply(counts(dds1), 1, gm_mean1)
dds1<-estimateSizeFactors(dds1, geoMeans = geoMeans1)

#Set reference level...............................................
dds1$Treatment<-relevel(dds1$Treatment,ref="Unburned")
levels(dds1$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dds1<-DESeq(dds1, fitType="local")


#Timepoint 2 **********************************************************************----
head(sample_data(physeqO3)$Treatment, 100) #look at the data

#Ran in case there are any n/a in data
physeqO3<-subset_samples(physeqO3, Treatment!= "None")

#Import phyloseq data to create a deseq object...................
dds2<-phyloseq_to_deseq2(physeqO3, ~ Treatment)

# calculate geometric means prior to estimate size factors.......
gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

geoMeans2 = apply(counts(dds2), 1, gm_mean2)
dds2<-estimateSizeFactors(dds2, geoMeans = geoMeans2)

#Set reference level.............................................
dds2$Treatment<-relevel(dds2$Treatment,ref="Unburned")
levels(dds2$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dds2<-DESeq(dds2, fitType="local")




#Timepoint 4 **********************************************************************----
head(sample_data(physeqO4)$Treatment, 100)

#Ran in case there are any n/a in data
physeqO4<-subset_samples(physeqO4, Treatment!= "None")

#Import phyloseq data to create a deseq object.....................
dds3<-phyloseq_to_deseq2(physeqO4, ~ Treatment)

# calculate geometric means prior to estimate size factors.........
gm_mean3 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans3 = apply(counts(dds3), 1, gm_mean3)
dds3<-estimateSizeFactors(dds3, geoMeans = geoMeans3)

#Set reference level..............................................
dds3$Treatment<-relevel(dds3$Treatment,ref="Unburned")
levels(dds3$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dds3<-DESeq(dds3, fitType="local")


#Timepoint 5 **********************************************************************----
head(sample_data(physeqO5)$Treatment, 100)

#Ran in case there are any n/a in data
physeqO5<-subset_samples(physeqO5, Treatment!= "None")

#Import phyloseq data to create a deseq object.......................
dds6<-phyloseq_to_deseq2(physeqO5, ~ Treatment)

# calculate geometric means prior to estimate size factors...........
gm_mean6 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans6 = apply(counts(dds6), 1, gm_mean6)
dds6<-estimateSizeFactors(dds6, geoMeans = geoMeans6)

#Set reference level.................................................
dds6$Treatment<-relevel(dds6$Treatment,ref="Unburned")
levels(dds6$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..................................
# (a.k.a. Gamma-Poisson) distribution
dds6<-DESeq(dds6, fitType="local")



#Timepoint 8 **********************************************************************----
head(sample_data(physeqO8)$Treatment, 100)

#Ran in case there are any n/a in data
physeqO8<-subset_samples(physeqO8, Treatment!= "None")

#Import phyloseq data to create a deseq object.......................
dds9<-phyloseq_to_deseq2(physeqO8, ~ Treatment)

# calculate geometric means prior to estimate size factors...........
gm_mean9 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans9 = apply(counts(dds9), 1, gm_mean9)
dds9<-estimateSizeFactors(dds9, geoMeans = geoMeans9)

#Set reference level.................................................
dds9$Treatment<-relevel(dds9$Treatment,ref="Unburned")
levels(dds9$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..................................
# (a.k.a. Gamma-Poisson) distribution
dds9<-DESeq(dds9, fitType="local")


#----
#----
#****************************************************************************************************-----
#--- RUN RESULT & EXPORT ------------------------------------------------------------------------------------------
#****************************************************************************************************-----
dir.create(("Fungi/Analysis/Deseq/Tables/"), recursive = T)

#Timepoint 1 ******************************************************************************************************************----
res1<-results(dds1)
res1<-res1[order(res1$padj, na.last=NA), ]
alpha<-0.05 #decrease this value and this should decrease the number of Genuss that you get 
sigtab1<-res1[(res1$padj < alpha), ];sigtab1
sigtab1<-cbind(as(sigtab1, "data.frame"), as(tax_table(physeqO2)[rownames(sigtab1), ], "matrix"));sigtab1
posigtab1<-sigtab1[sigtab1[, "log2FoldChange"] > 0, ]##Create ta table of most significant based on alpha and log2fold change.....

#Export results
write.csv(sigtab1,"Fungi/Analysis/Deseq/Tables/T1-Sigtab1-Original-Treatment-0.05.csv")
write.csv(posigtab1,"Fungi/Analysis/Deseq/Tables/T1-Posigtab1-Original-Treatment-0.05.csv")



#Timepoint 2 ***************************************************************************************************************----
res2<-results(dds2)
res2<-res2[order(res2$padj, na.last=NA), ]
alpha<-0.05
sigtab2<-res2[(res2$padj < alpha), ];sigtab2
sigtab2<-cbind(as(sigtab2, "data.frame"), as(tax_table(physeqO3)[rownames(sigtab2), ], "matrix"));sigtab2
posigtab2<-sigtab2[sigtab2[, "log2FoldChange"] > 0, ]


#Export results
write.csv(sigtab2,"Fungi/Analysis/Deseq/Tables/T2-Sigtab2-Original-Treatment-0.05.csv")
write.csv(posigtab2,"Fungi/Analysis/Deseq/Tables/T2-Posigtab2-Original-Treatment-0.05.csv")


#Timepoint 3 *******************************************************************************************----
res3<-results(dds3)
res3<-res3[order(res3$padj, na.last=NA), ]
alpha<-0.05
sigtab3<-res3[(res3$padj < alpha), ];sigtab3
sigtab3<-cbind(as(sigtab3, "data.frame"), as(tax_table(physeqO4)[rownames(sigtab3), ], "matrix"));sigtab3
posigtab3<-sigtab3[sigtab3[, "log2FoldChange"] > 0, ]

write.csv(sigtab3,"Fungi/Analysis/Deseq/Tables/T3-Sigtab3-Original-Treatment-0.05.csv")
write.csv(posigtab3,"Fungi/Analysis/Deseq/Tables/T3-Posigtab3-Original-Treatment-0.05.csv")

#Timepoint 6 *******************************************************************************************----
res6<-results(dds6)
res6<-res6[order(res6$padj, na.last=NA), ]
alpha<-0.05
sigtab6<-res6[(res6$padj < alpha), ];sigtab6
sigtab6<-cbind(as(sigtab6, "data.frame"), as(tax_table(physeqO5)[rownames(sigtab6), ], "matrix"));sigtab6
posigtab6<-sigtab6[sigtab6[, "log2FoldChange"] > 0, ]

write.csv(sigtab6,"Fungi/Analysis/Deseq/Tables/T6-Sigtab6-Original-Treatment-0.05.csv")
write.csv(posigtab6,"Fungi/Analysis/Deseq/Tables/T6-Posigtab6-Original-Treatment-0.05.csv")


#Timepoint 9 *********************************************************************************************----
res9<-results(dds9)
res9<-res9[order(res9$padj, na.last=NA), ]
alpha<-0.05
sigtab9<-res9[(res9$padj < alpha), ];sigtab9
sigtab9<-cbind(as(sigtab9, "data.frame"), as(tax_table(physeqO8)[rownames(sigtab9), ], "matrix"));sigtab9
posigtab9<-sigtab9[sigtab9[, "log2FoldChange"] > 0, ]

write.csv(sigtab9,"Fungi/Analysis/Deseq/Tables/T9-Sigtab9-Original-Treatment-0.05.csv")
write.csv(posigtab9,"Fungi/Analysis/Deseq/Tables/T9-Posigtab9-Original-Treatment-0.05.csv")


#----
#----
# Thawed samples.................................................................................................-----
#Timepoint 1 **********************************************************************----
head(sample_data(physeqT2)$Treatment, 100)

#Ran in case there are any n/a in data
physeqT2<-subset_samples(physeqT2, Treatment!= "None")

#Import phyloseq data to create a deseq object.....................
dds1<-phyloseq_to_deseq2(physeqT2, ~ Treatment)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans1 = apply(counts(dds1), 1, gm_mean1)
dds1<-estimateSizeFactors(dds1, geoMeans = geoMeans1)

#Set reference level...............................................
dds1$Treatment<-relevel(dds1$Treatment,ref="Unburned")
levels(dds1$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dds1<-DESeq(dds1, fitType="local")


#Timepoint 2 **********************************************************************----
head(sample_data(physeqT3)$Treatment, 100) #look at the data

#Ran in case there are any n/a in data
physeqT3<-subset_samples(physeqT3, Treatment!= "None")

#Import phyloseq data to create a deseq object...................
dds2<-phyloseq_to_deseq2(physeqT3, ~ Treatment)

# calculate geometric means prior to estimate size factors.......
gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

geoMeans2 = apply(counts(dds2), 1, gm_mean2)
dds2<-estimateSizeFactors(dds2, geoMeans = geoMeans2)

#Set reference level.............................................
dds2$Treatment<-relevel(dds2$Treatment,ref="Unburned")
levels(dds2$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dds2<-DESeq(dds2, fitType="local")




#Timepoint 4 **********************************************************************----
head(sample_data(physeqT4)$Treatment, 100)

#Ran in case there are any n/a in data
physeqT4<-subset_samples(physeqT4, Treatment!= "None")

#Import phyloseq data to create a deseq object.....................
dds3<-phyloseq_to_deseq2(physeqT4, ~ Treatment)

# calculate geometric means prior to estimate size factors.........
gm_mean3 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans3 = apply(counts(dds3), 1, gm_mean3)
dds3<-estimateSizeFactors(dds3, geoMeans = geoMeans3)

#Set reference level..............................................
dds3$Treatment<-relevel(dds3$Treatment,ref="Unburned")
levels(dds3$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dds3<-DESeq(dds3, fitType="local")


#Timepoint 5 **********************************************************************----
head(sample_data(physeqT5)$Treatment, 100)

#Ran in case there are any n/a in data
physeqT5<-subset_samples(physeqT5, Treatment!= "None")

#Import phyloseq data to create a deseq object.......................
dds6<-phyloseq_to_deseq2(physeqT5, ~ Treatment)

# calculate geometric means prior to estimate size factors...........
gm_mean6 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans6 = apply(counts(dds6), 1, gm_mean6)
dds6<-estimateSizeFactors(dds6, geoMeans = geoMeans6)

#Set reference level.................................................
dds6$Treatment<-relevel(dds6$Treatment,ref="Unburned")
levels(dds6$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..................................
# (a.k.a. Gamma-Poisson) distribution
dds6<-DESeq(dds6, fitType="local")



#Timepoint 8 **********************************************************************----
head(sample_data(physeqT8)$Treatment, 100)

#Ran in case there are any n/a in data
physeqT8<-subset_samples(physeqT8, Treatment!= "None")

#Import phyloseq data to create a deseq object.......................
dds9<-phyloseq_to_deseq2(physeqT8, ~ Treatment)

# calculate geometric means prior to estimate size factors...........
gm_mean9 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans9 = apply(counts(dds9), 1, gm_mean9)
dds9<-estimateSizeFactors(dds9, geoMeans = geoMeans9)

#Set reference level.................................................
dds9$Treatment<-relevel(dds9$Treatment,ref="Unburned")
levels(dds9$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..................................
# (a.k.a. Gamma-Poisson) distribution
dds9<-DESeq(dds9, fitType="local")


#--- RUN RESULT & EXPORT .......................................................................----
#Timepoint 1 ..........----
res1<-results(dds1)
res1<-res1[order(res1$padj, na.last=NA), ]
alpha<-0.05 #decrease this value and this should decrease the number of Genuss that you get 
sigtab1<-res1[(res1$padj < alpha), ];sigtab1
sigtab1<-cbind(as(sigtab1, "data.frame"), as(tax_table(physeqT2)[rownames(sigtab1), ], "matrix"));sigtab1
posigtab1<-sigtab1[sigtab1[, "log2FoldChange"] > 0, ]##Create ta table of most significant based on alpha and log2fold change.....

#Export results
write.csv(sigtab1,"Fungi/Analysis/Deseq/Tables/T1-Sigtab1-Thawed-Treatment-0.05.csv")
write.csv(posigtab1,"Fungi/Analysis/Deseq/Tables/T1-Posigtab1-Thawed-Treatment-0.05.csv")



#Timepoint 2 ***************************************************************************************************************----
res2<-results(dds2)
res2<-res2[order(res2$padj, na.last=NA), ]
alpha<-0.05
sigtab2<-res2[(res2$padj < alpha), ];sigtab2
sigtab2<-cbind(as(sigtab2, "data.frame"), as(tax_table(physeqT3)[rownames(sigtab2), ], "matrix"));sigtab2
posigtab2<-sigtab2[sigtab2[, "log2FoldChange"] > 0, ]


#Export results
write.csv(sigtab2,"Fungi/Analysis/Deseq/Tables/T2-Sigtab2-Thawed-Treatment-0.05.csv")
write.csv(posigtab2,"Fungi/Analysis/Deseq/Tables/T2-Posigtab2-Thawed-Treatment-0.05.csv")


#Timepoint 3 *******************************************************************************************----
res3<-results(dds3)
res3<-res3[order(res3$padj, na.last=NA), ]
alpha<-0.05
sigtab3<-res3[(res3$padj < alpha), ];sigtab3
sigtab3<-cbind(as(sigtab3, "data.frame"), as(tax_table(physeqT4)[rownames(sigtab3), ], "matrix"));sigtab3
posigtab3<-sigtab3[sigtab3[, "log2FoldChange"] > 0, ]

write.csv(sigtab3,"Fungi/Analysis/Deseq/Tables/T3-Sigtab3-Thawed-Treatment-0.05.csv")
write.csv(posigtab3,"Fungi/Analysis/Deseq/Tables/T3-Posigtab3-Thawed-Treatment-0.05.csv")

#Timepoint 6 *******************************************************************************************----
res6<-results(dds6)
res6<-res6[order(res6$padj, na.last=NA), ]
alpha<-0.05
sigtab6<-res6[(res6$padj < alpha), ];sigtab6
sigtab6<-cbind(as(sigtab6, "data.frame"), as(tax_table(physeqT5)[rownames(sigtab6), ], "matrix"));sigtab6
posigtab6<-sigtab6[sigtab6[, "log2FoldChange"] > 0, ]

write.csv(sigtab6,"Fungi/Analysis/Deseq/Tables/T6-Sigtab6-Thawed-Treatment-0.05.csv")
write.csv(posigtab6,"Fungi/Analysis/Deseq/Tables/T6-Posigtab6-Thawed-Treatment-0.05.csv")


#Timepoint 9 *********************************************************************************************----
res9<-results(dds9)
res9<-res9[order(res9$padj, na.last=NA), ]
alpha<-0.05
sigtab9<-res9[(res9$padj < alpha), ];sigtab9
sigtab9<-cbind(as(sigtab9, "data.frame"), as(tax_table(physeqT8)[rownames(sigtab9), ], "matrix"));sigtab9
posigtab9<-sigtab9[sigtab9[, "log2FoldChange"] > 0, ]

write.csv(sigtab9,"Fungi/Analysis/Deseq/Tables/T9-Sigtab9-Thawed-Treatment-0.05.csv")
write.csv(posigtab9,"Fungi/Analysis/Deseq/Tables/T9-Posigtab9-Thawed-Treatment-0.05.csv")









#---

#-----
#-----

# Frozen samples.................................................................................................-----
#Timepoint 1 **********************************************************************----
head(sample_data(physeqF2)$Treatment, 100)

#Ran in case there are any n/a in data
physeqF2<-subset_samples(physeqF2, Treatment!= "None")

#Import phyloseq data to create a deseq object.....................
dds1<-phyloseq_to_deseq2(physeqF2, ~ Treatment)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans1 = apply(counts(dds1), 1, gm_mean1)
dds1<-estimateSizeFactors(dds1, geoMeans = geoMeans1)

#Set reference level...............................................
dds1$Treatment<-relevel(dds1$Treatment,ref="Unburned")
levels(dds1$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dds1<-DESeq(dds1, fitType="local")


#Timepoint 2 **********************************************************************----
head(sample_data(physeqF3)$Treatment, 100) #look at the data

#Ran in case there are any n/a in data
physeqF3<-subset_samples(physeqF3, Treatment!= "None")

#Import phyloseq data to create a deseq object...................
dds2<-phyloseq_to_deseq2(physeqF3, ~ Treatment)

# calculate geometric means prior to estimate size factors.......
gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

geoMeans2 = apply(counts(dds2), 1, gm_mean2)
dds2<-estimateSizeFactors(dds2, geoMeans = geoMeans2)

#Set reference level.............................................
dds2$Treatment<-relevel(dds2$Treatment,ref="Unburned")
levels(dds2$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dds2<-DESeq(dds2, fitType="local")




#Timepoint 4 **********************************************************************----
head(sample_data(physeqF4)$Treatment, 100)

#Ran in case there are any n/a in data
physeqF4<-subset_samples(physeqF4, Treatment!= "None")

#Import phyloseq data to create a deseq object.....................
dds3<-phyloseq_to_deseq2(physeqF4, ~ Treatment)

# calculate geometric means prior to estimate size factors.........
gm_mean3 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans3 = apply(counts(dds3), 1, gm_mean3)
dds3<-estimateSizeFactors(dds3, geoMeans = geoMeans3)

#Set reference level..............................................
dds3$Treatment<-relevel(dds3$Treatment,ref="Unburned")
levels(dds3$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dds3<-DESeq(dds3, fitType="local")


#Timepoint 5 **********************************************************************----
head(sample_data(physeqF5)$Treatment, 100)

#Ran in case there are any n/a in data
physeqF5<-subset_samples(physeqF5, Treatment!= "None")

#Import phyloseq data to create a deseq object.......................
dds6<-phyloseq_to_deseq2(physeqF5, ~ Treatment)

# calculate geometric means prior to estimate size factors...........
gm_mean6 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans6 = apply(counts(dds6), 1, gm_mean6)
dds6<-estimateSizeFactors(dds6, geoMeans = geoMeans6)

#Set reference level.................................................
dds6$Treatment<-relevel(dds6$Treatment,ref="Unburned")
levels(dds6$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..................................
# (a.k.a. Gamma-Poisson) distribution
dds6<-DESeq(dds6, fitType="local")



#Timepoint 8 **********************************************************************----
head(sample_data(physeqF8)$Treatment, 100)

#Ran in case there are any n/a in data
physeqF8<-subset_samples(physeqF8, Treatment!= "None")

#Import phyloseq data to create a deseq object.......................
dds9<-phyloseq_to_deseq2(physeqF8, ~ Treatment)

# calculate geometric means prior to estimate size factors...........
gm_mean9 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans9 = apply(counts(dds9), 1, gm_mean9)
dds9<-estimateSizeFactors(dds9, geoMeans = geoMeans9)

#Set reference level.................................................
dds9$Treatment<-relevel(dds9$Treatment,ref="Unburned")
levels(dds9$Treatment)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..................................
# (a.k.a. Gamma-Poisson) distribution
dds9<-DESeq(dds9, fitType="local")


#--- RUN RESULT & EXPORT .......................................................................----
#Timepoint 1 ..........----
res1<-results(dds1)
res1<-res1[order(res1$padj, na.last=NA), ]
alpha<-0.05 #decrease this value and this should decrease the number of Genuss that you get 
sigtab1<-res1[(res1$padj < alpha), ];sigtab1
sigtab1<-cbind(as(sigtab1, "data.frame"), as(tax_table(physeqF2)[rownames(sigtab1), ], "matrix"));sigtab1
posigtab1<-sigtab1[sigtab1[, "log2FoldChange"] > 0, ]##Create ta table of most significant based on alpha and log2fold change.....

#Export results
write.csv(sigtab1,"Fungi/Analysis/Deseq/Tables/T1-Sigtab1-Frozen-Treatment-0.05.csv")
write.csv(posigtab1,"Fungi/Analysis/Deseq/Tables/T1-Posigtab1-Frozen-Treatment-0.05.csv")



#Timepoint 2 ***************************************************************************************************************----
res2<-results(dds2)
res2<-res2[order(res2$padj, na.last=NA), ]
alpha<-0.05
sigtab2<-res2[(res2$padj < alpha), ];sigtab2
sigtab2<-cbind(as(sigtab2, "data.frame"), as(tax_table(physeqF3)[rownames(sigtab2), ], "matrix"));sigtab2
posigtab2<-sigtab2[sigtab2[, "log2FoldChange"] > 0, ]


#Export results
write.csv(sigtab2,"Fungi/Analysis/Deseq/Tables/T2-Sigtab2-Frozen-Treatment-0.05.csv")
write.csv(posigtab2,"Fungi/Analysis/Deseq/Tables/T2-Posigtab2-Frozen-Treatment-0.05.csv")


#Timepoint 3 *******************************************************************************************----
res3<-results(dds3)
res3<-res3[order(res3$padj, na.last=NA), ]
alpha<-0.05
sigtab3<-res3[(res3$padj < alpha), ];sigtab3
sigtab3<-cbind(as(sigtab3, "data.frame"), as(tax_table(physeqF4)[rownames(sigtab3), ], "matrix"));sigtab3
posigtab3<-sigtab3[sigtab3[, "log2FoldChange"] > 0, ]

write.csv(sigtab3,"Fungi/Analysis/Deseq/Tables/T3-Sigtab3-Frozen-Treatment-0.05.csv")
write.csv(posigtab3,"Fungi/Analysis/Deseq/Tables/T3-Posigtab3-Frozen-Treatment-0.05.csv")

#Timepoint 6 *******************************************************************************************----
res6<-results(dds6)
res6<-res6[order(res6$padj, na.last=NA), ]
alpha<-0.05
sigtab6<-res6[(res6$padj < alpha), ];sigtab6
sigtab6<-cbind(as(sigtab6, "data.frame"), as(tax_table(physeqF5)[rownames(sigtab6), ], "matrix"));sigtab6
posigtab6<-sigtab6[sigtab6[, "log2FoldChange"] > 0, ]

write.csv(sigtab6,"Fungi/Analysis/Deseq/Tables/T6-Sigtab6-Frozen-Treatment-0.05.csv")
write.csv(posigtab6,"Fungi/Analysis/Deseq/Tables/T6-Posigtab6-Frozen-Treatment-0.05.csv")


#Timepoint 9 *********************************************************************************************----
res9<-results(dds9)
res9<-res9[order(res9$padj, na.last=NA), ]
alpha<-0.05
sigtab9<-res9[(res9$padj < alpha), ];sigtab9
sigtab9<-cbind(as(sigtab9, "data.frame"), as(tax_table(physeqF8)[rownames(sigtab9), ], "matrix"));sigtab9
posigtab9<-sigtab9[sigtab9[, "log2FoldChange"] > 0, ]

write.csv(sigtab9,"Fungi/Analysis/Deseq/Tables/T9-Sigtab9-Frozen-Treatment-0.05.csv")
write.csv(posigtab9,"Fungi/Analysis/Deseq/Tables/T9-Posigtab9-Frozen-Treatment-0.05.csv")









#---





#----
#****************************************************************************************************-----
#---CREATE GRAPHS---------------------------------------------------------------------------------------
#****************************************************************************************************-----
#Check unique names
#Reimport files...........................................................................................----
CountsAllO<-read.csv("Fungi/Analysis/Deseq/Tables/Deseq-All-Counts-Original.csv")
CountsAllT<-read.csv("Fungi/Analysis/Deseq/Tables/All-Thawed-Barplot.csv")
CountsAllF<-read.csv("Fungi/Analysis/Deseq/Tables/All-Frozen-Barplot.csv")


#Create a graph pf all the timepoints..........................................................-------
attach(CountsAll0)

Vertical<-ggplot(data=CountsAll0,  aes(x=reorder(Genus, +Counts2), y=Counts2, fill=as.factor(TSF)))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = c("#a6cce0","#f0ab43","#8a78ad","#286a9c","#630f13")) +
  #scale_fill_manual(values = c("#1c5a75","#758c96","#d3e1e6","#dbcba7","#876354"))+
  labs(title = "(a) Fungi Original", tag = "", face="bold") +
  geom_hline(yintercept=0, color = "#4b4c4d", linewidth=.9, linetype = "dotted")+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=12,face="bold"),
        axis.title = element_text(size = 16, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.text.x=element_text(size=14, colour = "black"),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=14),
        legend.position = "right")+
  labs(x = "", y="Number of Differentially Abundant Genuss")+
  coord_flip(ylim = c(-5,5))+
  scale_y_continuous(limits = c(-5,5), 
                     breaks = c(-5,-2.5,-0, 2.5, 5))+
  guides(fill=guide_legend(ncol =1))+
  labs(fill = "Functions");Vertical


pdf("Fungi/Analysis/Deseq/All-Original-Barplot.pdf", height=4, width=9)
Vertical
dev.off()

#Thawed samples
VerticalT<-ggplot(data=CountsAllT,  aes(x=reorder(Genus, +Counts2), y=Counts2, fill=as.factor(TSF)))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = c("#a6cce0","#f0ab43","#8a78ad","#286a9c","#630f13")) +
  #scale_fill_manual(values = c("#1c5a75","#758c96","#d3e1e6","#dbcba7","#876354"))+
  labs(title = "(a) Fungi Thawed", tag = "", face="bold") +
  geom_hline(yintercept=0, color = "#4b4c4d", linewidth=.9, linetype = "dotted")+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=12,face="bold"),
        axis.title = element_text(size = 16, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.text.x=element_text(size=14, colour = "black"),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=14),
        legend.position = "right")+
  labs(x = "", y="Number of Differentially Abundant Genuss")+
  coord_flip(ylim = c(-5,5))+
  scale_y_continuous(limits = c(-5,5), 
                     breaks = c(-5,-2.5,-0, 2.5, 5))+
  guides(fill=guide_legend(ncol =1))+
  labs(fill = "Functions");VerticalT


pdf("Fungi/Analysis/Deseq/All-Thawed-Barplot.pdf", height=4, width=9)
VerticalT
dev.off()


#Frozen samples
VerticalF<-ggplot(data=CountsAllF,  aes(x=reorder(Genus, +Counts2), y=Counts2, fill=as.factor(TSF)))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = c("#a6cce0","#f0ab43","#8a78ad","#286a9c","#630f13")) +
  #scale_fill_manual(values = c("#1c5a75","#758c96","#d3e1e6","#dbcba7","#876354"))+
  labs(title = "(a) Fungi Frozen", tag = "", face="bold") +
  geom_hline(yintercept=0, color = "#4b4c4d", linewidth=.9, linetype = "dotted")+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=12,face="bold"),
        axis.title = element_text(size = 16, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.text.x=element_text(size=14, colour = "black"),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=14),
        legend.position = "right")+
  labs(x = "", y="Number of Differentially Abundant Genuss")+
  coord_flip(ylim = c(-5,5))+
  scale_y_continuous(limits = c(-5,5), 
                     breaks = c(-5,-2.5,-0, 2.5, 5))+
  guides(fill=guide_legend(ncol =1))+
  labs(fill = "Functions");VerticalF


pdf("Fungi/Analysis/Deseq/All-Frozen-Barplot.pdf", height=4, width=9)
VerticalF
dev.off()



CorAll<-Vertical / VerticalT / VerticalF +
  plot_layout(ncol = 3, nrow = 1) &
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom');CorAll


pdf("Fungi/Analysis/Deseq/Deseq-Fungi-All.pdf", height=11, width=14)
CorAll
dev.off()
