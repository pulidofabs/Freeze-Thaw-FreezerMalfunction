#Reset R's Brain
rm(list=ls())

#Set working directory................................................................................................
setwd("~/Dropbox/2-UC-Riverside/1-Dissertation/1-Research/1-HolyFire/4-Freeze-thaw-Analysis/1-Lib7Thawed-NoCNF9WT2/")

#Load Packages-------------------------------------------------------------------------------
library("qiime2R")
library(DESeq2)
library(ggplot2)#Plotting
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(ochRe)# color paletteibrary(multcomp)#multiple comparisons--for glm 
library(tidyverse)
library(tidyr)

Metadata<-read_tsv("Bacteria/Metadata/MetaRareBac-CNF9WT2.tsv")#Rare metadata
Table<-read_qza("Qiime/Merged/Bacteria/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")
Tree<-read_qza("Qiime/Merged/Bacteria/CoreMetrics/Phylo/Bacteria-rooted-tree-1-4RT-9WT2.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Merged/Bacteria/Taxonomy/Bacteria-taxonomy-1-4RT-9WT2-132.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                                                     c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#Create phyloseq artifact.........................................................................
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleIDLib")))





#----
#----
##----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("_D0__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("_D1__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("_D2__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("_D3__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("_D4__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("_D5__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("_D6__", "", tax_table(physeq)[, "Species"])


sample_data(physeq)$SampleType<-factor(sample_data(physeq)$SampleType) # make as factor

#Subset data to look at sample types independently ----------------------------------------------------------------------
physeqOT<-subset_samples(physeq, SampleType != "FrozenDNA");sample_data(physeqOT)$SampleType
physeqOF<-subset_samples(physeq, SampleType != "Thawed");sample_data(physeqOF)$SampleType
physeqTF<-subset_samples(physeq, SampleType != "Original");sample_data(physeqTF)$SampleType



#----
#-----
#****************************************************************************************************-----
#--- RUN DESEQ FOR ORIGINAL SAMPLES
#****************************************************************************************************-----

#Original vs Thawed samples.......................................................----
head(sample_data(physeqOT)$SampleType, 100)

#Ran in case there are any n/a in data
physeqOT<-subset_samples(physeqOT, SampleType!= "None")

#Import phyloseq data to create a deseq object.....................
ddsOT<-phyloseq_to_deseq2(physeqOT, ~ SampleType)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans1 = apply(counts(ddsOT), 1, gm_mean1)
ddsOT<-estimateSizeFactors(ddsOT, geoMeans = geoMeans1)

#Set reference level...............................................
ddsOT$SampleType<-relevel(ddsOT$SampleType,ref="Original")
levels(ddsOT$SampleType)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
ddsOT<-DESeq(ddsOT, fitType="local")
#-----

#Original vs Frozen samples.......................................................----
head(sample_data(physeqOF)$SampleType, 100)

#Ran in case there are any n/a in data
physeqOF<-subset_samples(physeqOF, SampleType!= "None")

#Import phyloseq data to create a deseq object.....................
ddsOF<-phyloseq_to_deseq2(physeqOF, ~ SampleType)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans1 = apply(counts(ddsOF), 1, gm_mean1)
ddsOF<-estimateSizeFactors(ddsOF, geoMeans = geoMeans1)

#Set reference level...............................................
ddsOF$SampleType<-relevel(ddsOF$SampleType,ref="Original")
levels(ddsOF$SampleType)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
ddsOF<-DESeq(ddsOF, fitType="local")

#-----

#Thawed vs Frozen samples.......................................................----
head(sample_data(physeqTF)$SampleType, 100)

#Ran in case there are any n/a in data
physeqTF<-subset_samples(physeqTF, SampleType!= "None")

#Import phyloseq data to create a deseq object.....................
ddsTF<-phyloseq_to_deseq2(physeqTF, ~ SampleType)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans1 = apply(counts(ddsTF), 1, gm_mean1)
ddsTF<-estimateSizeFactors(ddsTF, geoMeans = geoMeans1)

#Set reference level...............................................
ddsTF$SampleType<-relevel(ddsTF$SampleType,ref="FrozenDNA")
levels(ddsTF$SampleType)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
ddsTF<-DESeq(ddsTF, fitType="local")






#----
#----
#****************************************************************************************************-----
#--- RUN RESULT & EXPORT ------------------------------------------------------------------------------------------
#****************************************************************************************************-----
dir.create(("Bacteria/Analysis/Deseq/Tables/"), recursive = T)

#Original-frozen........................................................................................-----
resOF<-results(ddsOF)
resOF<-resOF[order(resOF$padj, na.last=NA), ]
alpha<-0.05 #decrease this value and this should decrease the number of Genuss that you get 
sigtabOF<-resOF[(resOF$padj < alpha), ];sigtabOF
sigtabOF<-cbind(as(sigtabOF, "data.frame"), as(tax_table(physeqOT)[rownames(sigtabOF), ], "matrix"));sigtabOF
posigtabOF<-sigtabOF[sigtabOF[, "log2FoldChange"] > 0, ]##Create ta table of most significant based on alpha and log2fold change.....

#Export results
write.csv(sigtabOF,"Bacteria/Analysis/Deseq/Tables/T1-Sigtab-Original-Frozen-0.05.csv")
write.csv(posigtabOF,"Bacteria/Analysis/Deseq/Tables/T1-Posigtab1-Original-Frozen-0.05.csv")


#Original-thawed........................................................................................-----
resOT<-results(ddsOT)
resOT<-resOT[order(resOT$padj, na.last=NA), ]
alpha<-0.05 #decrease this value and this should decrease the number OT Genuss that you get 
sigtabOT<-resOT[(resOT$padj < alpha), ];sigtabOT
sigtabOT<-cbind(as(sigtabOT, "data.frame"), as(tax_table(physeqOT)[rownames(sigtabOT), ], "matrix"));sigtabOT
posigtabOT<-sigtabOT[sigtabOT[, "log2FoldChange"] > 0, ]##Create ta table OT most significant based on alpha and log2fold change.....

#Export results
write.csv(sigtabOT,"Bacteria/Analysis/Deseq/Tables/T1-Sigtab-Original-Thawed-0.05.csv")
write.csv(posigtabOT,"Bacteria/Analysis/Deseq/Tables/T1-Posigtab1-Original-Thawed-0.05.csv")


#Frozen-thawed........................................................................................-----
resTF<-results(ddsTF)
resTF<-resTF[order(resTF$padj, na.last=NA), ]
alpha<-0.05 #decrease this value and this should decrease the number TF Genuss that you get 
sigtabTF<-resTF[(resTF$padj < alpha), ];sigtabTF
sigtabTF<-cbind(as(sigtabTF, "data.frame"), as(tax_table(physeqTF)[rownames(sigtabTF), ], "matrix"));sigtabTF
posigtabTF<-sigtabTF[sigtabTF[, "log2FoldChange"] > 0, ]##Create ta table TF most significant based on alpha and log2fold change.....

#Export results
write.csv(sigtabTF,"Bacteria/Analysis/Deseq/Tables/T1-Sigtab-Frozen-Thawed-0.05.csv")
write.csv(posigtabTF,"Bacteria/Analysis/Deseq/Tables/T1-Posigtab1-Frozen-Thawed-0.05.csv")


#No significance in any case

#----
#****************************************************************************************************-----
#---CREATE GRAPHS---------------------------------------------------------------------------------------
#****************************************************************************************************-----
#Check unique names
#Reimport files...........................................................................................----
CountsAllO<-read.csv("Bacteria/Analysis/Deseq/Tables/All-Original-Barplot.csv")
CountsAllT<-read.csv("Bacteria/Analysis/Deseq/Tables/All-Thawed-Barplot.csv")
CountsAllF<-read.csv("Bacteria/Analysis/Deseq/Tables/All-Frozen-Barplot.csv")


#Create a graph pf all the timepoints..........................................................-------
attach(CountsAll0)

Vertical<-ggplot(data=CountsAllO,  aes(x=reorder(Genus, +Counts2), y=Counts2, fill=as.factor(TSF)))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = c("#a6cce0","#f0ab43","#8a78ad","#286a9c","#630f13")) +
  #scale_fill_manual(values = c("#1c5a75","#758c96","#d3e1e6","#dbcba7","#876354"))+
  labs(title = "(a) Bacteria Original", tag = "", face="bold") +
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
  labs(x = "", y="Number of Differentially Abundant Genuss") +
  coord_flip(ylim = c(-5,13))+
  scale_y_continuous(limits = c(-5,13), 
                     breaks = c(-5,-2.5,-0, 2.5, 5,10,12.5))+
  guides(fill=guide_legend(ncol =1))+
  labs(fill = "Bacctions");Vertical


pdf("Bacteria/Analysis/Deseq/All-Original-Barplot.pdf", height=4, width=9)
Vertical
dev.off()

#Thawed samples
VerticalT<-ggplot(data=CountsAllT,  aes(x=reorder(Genus, +Counts2), y=Counts2, fill=as.factor(TSF)))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = c("#a6cce0","#f0ab43","#8a78ad","#286a9c","#630f13")) +
  #scale_fill_manual(values = c("#1c5a75","#758c96","#d3e1e6","#dbcba7","#876354"))+
  labs(title = "(a) Bacteria Thawed", tag = "", face="bold") +
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
  labs(x = "", y="Number of Differentially Abundant Genuss") +
  coord_flip(ylim = c(-10,10))+
  scale_y_continuous(limits = c(-10,10), 
                     breaks = c(-10,-5,-0, 5, 10))+
  guides(fill=guide_legend(ncol =1))+
  labs(fill = "Bacctions");VerticalT


pdf("Bacteria/Analysis/Deseq/All-Thawed-Barplot.pdf", height=4, width=9)
VerticalT
dev.off()


#Frozen samples
VerticalF<-ggplot(data=CountsAllF,  aes(x=reorder(Genus, +Counts2), y=Counts2, fill=as.factor(TSF)))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = c("#a6cce0","#f0ab43","#8a78ad","#286a9c","#630f13")) +
  #scale_fill_manual(values = c("#1c5a75","#758c96","#d3e1e6","#dbcba7","#876354"))+
  labs(title = "(a) Bacteria Frozen", tag = "", face="bold") +
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
  coord_flip(ylim = c(-10,10))+
  scale_y_continuous(limits = c(-10,10), 
                     breaks = c(-10,-5,-0, 5, 10))+
  guides(fill=guide_legend(ncol =1))+
  labs(fill = "Bacctions");VerticalF


pdf("Bacteria/Analysis/Deseq/All-Frozen-Barplot.pdf", height=4, width=9)
VerticalF
dev.off()



CorAll<-Vertical / VerticalT / VerticalF +
  plot_layout(ncol = 3, nrow = 1) &
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom');CorAll


pdf("Bacteria/Analysis/Deseq/Deseq-Bacteria-All.pdf", height=11, width=25)
CorAll
dev.off()
