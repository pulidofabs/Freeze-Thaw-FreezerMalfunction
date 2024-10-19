#Date last ran
# 5/12/2022

#Reset R's 
rm(list=ls())

#Set working directory................................................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/Lib7Thawed-NoCNF9WT2/")


#Load librarires......................................................................................................
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)#to load qiime object
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(patchwork)#create graph panels
library(ochRe)
library(patchwork)

#----
#----
#*********************************************************************************************************----
# DATA FOR PHYLOSEQ AND CALC ABUNDANCE (ONLY DO THE 1ST TIME, AFTER JUST LOAD DATA ON LINE 190 ) -------------------------------------------------------
#*********************************************************************************************************----
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
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----

rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("D_0__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("D_1__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("D_2__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("D_3__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("D_4__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("D_5__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("D_6__", "", tax_table(physeq)[, "Species"])


#--Subset data by treatment (burned vs unburned) ...................................................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor

# * * Subset All to maintain only the burned samples of all data...................................................
physeqB<-subset_samples(physeq,Treatment=="Burned");physeqB#7926x30
physeqUn<-subset_samples(physeq,Treatment=="Unburned");physeqUn#7926x24

#Subset burned samples by sample type so that we can see if we can see TSF...............
PhyThawB<-subset_samples(physeqB,SampleType=="Thawed");PhyThawB#7926x10
PhyFrozB<-subset_samples(physeqB,SampleType=="FrozenDNA");PhyFrozB#7926x10
PhyOrigB<-subset_samples(physeqB,SampleType=="Original");PhyOrigB#7926x10

#Subset Unburned samples by sample type so that we can see if we can see TSF...............
PhyThawUn<-subset_samples(physeqUn,SampleType=="Thawed");PhyThawUn#7926x9
PhyFrozUn<-subset_samples(physeqUn,SampleType=="FrozenDNA");PhyFrozUn#7926x7
PhyOrigUn<-subset_samples(physeqUn,SampleType=="Original");PhyOrigUn#7926x8




#----
#-----
#************************************************************************************************************------
#*----------------------CALCULATE RELATIVE ABUDNANCE -------------------------------------------------------------
#************************************************************************************************************------
#** Sample Tupe analysis Burned samples............................................................................
GenBSt <- tax_glom(physeqB, taxrank = 'Genus')
(GenBSt= merge_samples(GenBSt , "SampleType")) #merge them by variable of intest
RelGenBSt <- transform_sample_counts(GenBSt, function(x) x/sum(x)) #calculate rel abund
RelGenBSt1 <- psmelt(RelGenBSt) #create a dataframe to work in R
RelGenBSt1$Genus <- as.character(RelGenBSt1$Genus)
RelGenBSt1$Genus[RelGenBSt1$Abundance < 0.03] <- "< 3% abund." 


#** Stple Type unburned analysis.................................................................................
GenUSt <- tax_glom(physeqUn, taxrank = 'Genus')
(GenUSt= merge_samples(GenUSt , "SampleType")) #merge them by variable of intest
RelGenUSt <- transform_sample_counts(GenUSt, function(x) x/sum(x)) #calculate rel abund
RelGenUSt1 <- psmelt(RelGenUSt) #create a dataframe to work in R
RelGenUSt1$Genus <- as.character(RelGenUSt1$Genus)
RelGenUSt1$Genus[RelGenUSt1$Abundance < 0.03] <- "< 3% abund." 


#TSF Burned analysis.................................................................................-------
#Look at each individual sample type by TSF to see if anything changes..................................-----
#Thawed samples........................................................
GenTB <- tax_glom(PhyThawB, taxrank = 'Genus')
(GenTBtsf  = merge_samples(GenTB, "TSFdays")) 
sample_data(GenTBtsf)$TSFdays <- factor(sample_names(GenTBtsf))
RelGenTBtsf <- transform_sample_counts(GenTBtsf, function(x) x/sum(x)) 
RelGenTBtsf1 <- psmelt(RelGenTBtsf)
RelGenTBtsf1$Genus <- as.character(RelGenTBtsf1$Genus)
RelGenTBtsf1$Genus[RelGenTBtsf1$Abundance < 0.03] <- "< 3% abund."

#Original samples.............................................................
GenOB <- tax_glom(PhyOrigB, taxrank = 'Genus')
(GenOBtsf  = merge_samples(GenOB, "TSFdays")) 
sample_data(GenOBtsf)$TSFdays <- factor(sample_names(GenOBtsf))
RelGenOBtsf <- transform_sample_counts(GenOBtsf, function(x) x/sum(x)) 
RelGenOBtsf1 <- psmelt(RelGenOBtsf)
RelGenOBtsf1$Genus <- as.character(RelGenOBtsf1$Genus)
RelGenOBtsf1$Genus[RelGenOBtsf1$Abundance < 0.03] <- "< 3% abund."

#Frozen DNA samples........................................................
GenFB <- tax_glom(PhyFrozB, taxrank = 'Genus')
(GenFBtsf  = merge_samples(GenFB, "TSFdays")) 
sample_data(GenFBtsf)$TSFdays <- factor(sample_names(GenFBtsf))
RelGenFBtsf <- transform_sample_counts(GenFBtsf, function(x) x/sum(x)) 
RelGenFBtsf1 <- psmelt(RelGenFBtsf)
RelGenFBtsf1$Genus <- as.character(RelGenFBtsf1$Genus)
RelGenFBtsf1$Genus[RelGenFBtsf1$Abundance < 0.03] <- "< 3% abund."


#TSF Unburned analysis.................................................................................------
#Look at each individual sample type by TSF to see if anything changes..................................
#Thawed samples............................................................
GenTU <- tax_glom(PhyThawUn, taxrank = 'Genus')
(GenTUtsf  = merge_samples(GenTU, "TSFdays")) 
sample_data(GenTUtsf)$TSFdays <- factor(sample_names(GenTUtsf))
RelGenTUtsf <- transform_sample_counts(GenTUtsf, function(x) x/sum(x)) 
RelGenTUtsf1 <- psmelt(RelGenTUtsf)
RelGenTUtsf1$Genus <- as.character(RelGenTUtsf1$Genus)
RelGenTUtsf1$Genus[RelGenTUtsf1$Abundance < 0.03] <- "< 3% abund."


#Original samples.........................................................
GenOU <- tax_glom(PhyOrigUn, taxrank = 'Genus')
(GenOUtsf  = merge_samples(GenOU, "TSFdays")) 
sample_data(GenOUtsf)$TSFdays <- factor(sample_names(GenOUtsf))
RelGenOUtsf <- transform_sample_counts(GenOUtsf, function(x) x/sum(x)) 
RelGenOUtsf1 <- psmelt(RelGenOUtsf)
RelGenOUtsf1$Genus <- as.character(RelGenOUtsf1$Genus)
RelGenOUtsf1$Genus[RelGenOUtsf1$Abundance < 0.03] <- "< 3% abund."

#Frozen DNA samples........................................................
GenFU <- tax_glom(PhyFrozUn, taxrank = 'Genus')
(GenFUtsf  = merge_samples(GenFU, "TSFdays")) 
sample_data(GenFUtsf)$TSFdays <- factor(sample_names(GenFUtsf))
RelGenFUtsf <- transform_sample_counts(GenFUtsf, function(x) x/sum(x)) 
RelGenFUtsf1 <- psmelt(RelGenFUtsf)
RelGenFUtsf1$Genus <- as.character(RelGenFUtsf1$Genus)
RelGenFUtsf1$Genus[RelGenFUtsf1$Abundance < 0.03] <- "< 3% abund."



#Export files so that we do not need to rerun above to fix the graphs...................................-----
dir.create("Bacteria/Analysis/RelativeAbundance/Tables/1Per/", recursive = TRUE)#create a directory
dir.create("Bacteria/Analysis/RelativeAbundance/Tables/3Per/", recursive = TRUE)#create a directory
dir.create("Bacteria/Analysis/RelativeAbundance/Tables/2Per/", recursive = TRUE)#create a directory

#3 percent data ....................................................................................................
write.csv(RelGenBSt1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/SampleType-RelAbund-3-9WT2-BurnBac.csv")
write.csv(RelGenUSt1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/SampleType-RelAbund-3-9WT2-UnBac.csv")

write.csv(RelGenTBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/ThawBurn-RelAbund-3-9WT2-Bac.csv")
write.csv(RelGenOBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/OrigBurn-RelAbund-3-9WT2-Bac.csv")
write.csv(RelGenFBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/FrozenBurn-RelAbund-3-9WT2-Bac.csv")

write.csv(RelGenTUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/ThawUnburn-RelAbund-3-9WT2-Bac.csv")
write.csv(RelGenOUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/OrigUnburn-RelAbund-3-9WT2-Bac.csv")
write.csv(RelGenFUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/3Per/FrozenUnburn-RelAbund-3-9WT2-Bac.csv")


#----
#------
#2 percent data ....................................................................................................
write.csv(RelGenBSt1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/SampleType-RelAbund-2-9WT2-BurnBac.csv")
write.csv(RelGenUSt1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/SampleType-RelAbund-2-9WT2-UnBac.csv")

write.csv(RelGenTBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/ThawBurn-RelAbund-2-9WT2-Bac.csv")
write.csv(RelGenOBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/OrigBurn-RelAbund-2-9WT2-Bac.csv")
write.csv(RelGenFBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/FrozenBurn-RelAbund-2-9WT2-Bac.csv")

write.csv(RelGenTUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/ThawUnburn-RelAbund-2-9WT2-Bac.csv")
write.csv(RelGenOUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/OrigUnburn-RelAbund-2-9WT2-Bac.csv")
write.csv(RelGenFUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/2Per/FrozenUnburn-RelAbund-2-9WT2-Bac.csv")



#----
#------
#1 percent data ....................................................................................................
write.csv(RelGenBSt1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/SampleType-RelAbund-1-9WT2-BurnBac.csv")
write.csv(RelGenUSt1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/SampleType-RelAbund-1-9WT2-UnBac.csv")

write.csv(RelGenTBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/ThawBurn-RelAbund-1-9WT2-Bac.csv")
write.csv(RelGenOBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/OrigBurn-RelAbund-1-9WT2-Bac.csv")
write.csv(RelGenFBtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/FrozenBurn-RelAbund-1-9WT2-Bac.csv")

write.csv(RelGenTUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/ThawUnburn-RelAbund-1-9WT2-Bac.csv")
write.csv(RelGenOUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/OrigUnburn-RelAbund-1-9WT2-Bac.csv")
write.csv(RelGenFUtsf1,"Bacteria/Analysis/RelativeAbundance/Tables/1Per/FrozenUnburn-RelAbund-1-9WT2-Bac.csv")



#Look at your exported files and if there are any uncultured samples, manually blast and 
#see if you can get a better name, if not, use the furthest taxonomy level assigned and 
#change your excel sheet to show this name (for plotting), better to have a name than
#"uncultured or unknown"



#----
#----
#**************************************************************************************************************----
# -lLOAD DATA AND QUALITY CONTROL FOR GRAPHING PURPOSES -------------------------------------------------
#**************************************************************************************************************----
#Reimport files to create graphs (save you time from having to run the above)
#3% level identification................................................................................
RelGenBSt3<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/3Per/SampleType-RelAbund-3-9WT2-BurnBac.csv")
RelGenUSt3<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/3Per/SampleType-RelAbund-3-9WT2-UnBac.csv")

#2% level identification.................................................................................
RelGenBSt2<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/2Per/SampleType-RelAbund-2-9WT2-BurnBac.csv")
RelGenUSt2<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/2Per/SampleType-RelAbund-2-9WT2-UnBac.csv")

#1% level identification.................................................................................
RelGenBSt1<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/1Per/SampleType-RelAbund-1-9WT2-BurnBac.csv")
RelGenUSt1<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/1Per/SampleType-RelAbund-1-9WT2-UnBac.csv")


#TSF-Graphs......----------------------------------------------------------------------------------------------------
#RelGenTBtsf1<-read.csv("Analysis/RelativeAbundance/Tables/ThawBurn-RelAbund-3-Fun.csv")
#RelGenOBtsf1<-read.csv("Analysis/RelativeAbundance/Tables/OrigBurn-RelAbund-3-Fun.csv")
#RelGenFBtsf1<-read.csv("Analysis/RelativeAbundance/Tables/FrozenBurn-RelAbund-3-Fun.csv")

#RelGenTUtsf1<-read.csv("Analysis/RelativeAbundance/Tables/ThawUnburn-RelAbund-3-Fun.csv")
#RelGenOUtsf1<-read.csv("Analysis/RelativeAbundance/Tables/OrigUnburn-RelAbund-3-Fun.csv")
#RelGenFUtsf1<-read.csv("Analysis/RelativeAbundance/Tables/FrozenUnburn-RelAbund-3-Fun.csv")


#-----
#Reorder Time Since fire .....................................................................................................
#Genus TSF Unburned samples
RelGenTUtsf1$Sample <- factor(RelGenTUtsf1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))
RelGenOUtsf1$Sample <- factor(RelGenOUtsf1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))
RelGenFUtsf1$Sample <- factor(RelGenFUtsf1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))

#Genus TSF Burned samples
RelGenTBtsf1$Sample <- factor(RelGenTBtsf1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))
RelGenOBtsf1$Sample <- factor(RelGenOBtsf1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))
RelGenFBtsf1$Sample <- factor(RelGenFBtsf1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))

#Look at length so that you know how many colors you need to create...............................................-----
#look at the names to see if you need to rename any of them..............................................
length(unique(RelGenBSt3$Genus));unique(RelGenBSt3$Genus)#(5)
length(unique(RelGenBSt1$Genus));unique(RelGenBSt1$Genus)#(19)
length(unique(RelGenBSt2$Genus));unique(RelGenBSt2$Genus)#(6)


#-
#-----
#Rename the bacteria taxa (name is too long for plotting)...........................................................
#Treatment
RelGenBSt3$Genus[RelGenBSt3$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 
RelGenUSt3$Genus[RelGenUSt3$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 
RelGenUSt3$Genus[RelGenUSt3$Genus =="Burkholderia-Caballeronia-Paraburkholderia" ]<-"Burkholderia" 



#
#Create color scheme for the graphs................................................................................
Trt<-c("< 3% abund."="#dee0df","Bryobacter"="#baa4bf",
         "Burkholderia"="#d4c9d6","Can. Udaeobacter"="#faf7d2",
         "Clostridiales"="#e0d4b1","Domibacillus"="#d6c07e",
         "Massilia"="#a3a27e","Mycobacterium"="#787612",
         "Paenibacillus"="#425778","RB41"="#96806f")



TrtUn<-c("< 3% abund."="#dee0df", "Acetobacteraceae" = "#f5d9d5",
         "Acidibacter"="#c49993","Acidipila"="#9e6057","Actinoplanes" = "#a11b1b",
         "Blastococcus"="#6a5170","Bradyrhizobium"="#9c7fa3",
         "Bryobacter"="#baa4bf","Burkholderia"="#d4c9d6","Can. Udaeobacter"="#faf7d2",
         "Candidatus Solibacter" ="#fff894", "Chitinophagaceae"="#c7c05f",
         "Chthoniobacter"="#827c2f", "Conexibacter"="#e6bc57",
         "Gemmatimonadaceae"="#c4a64b","Gemmatimonas"="#d1d0b4",
         "Granulicella"= "#acbda4","Haliangium"="#5a944d","Kribbella" = "#267813",
         "Massilia"="#145405","Micropepsaceae"="#a3a27e", "Mucilaginibacter"="#63633b",
         "Mycobacterium"="#787612","Pseudonocardia"="#b0a094", "Phenylobacterium"="#c7d4c5",
         "RB41"="#96806f", "Solirubrobacter"="#664e3c","Sphingomonas"="#d3e2e3",
         "Streptomyces"="#272829","Tepidisphaerales"="#93b8ba","Xanthobacteraceae" ="#5a9396")









#----
#----
#*********************************************************************************************************----
# --------------------------------- TREATMENT PER SAMPLE TYPE (1%)-------------------------------------------------------
#*********************************************************************************************************----
SampB<-ggplot(data=RelGenBSt3, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack",show.legend = T) + #,show.legend = F
  #scale_fill_ochre(palette="olsen_seq")
  scale_fill_manual(values = Trt)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        #axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
        legend.position="right",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=20))+ 
  guides(fill=guide_legend(ncol = 1))+
  xlab("")+ylab("");SampB

SampUn<-ggplot(data=RelGenUSt3, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack", show.legend = F) + 
  #scale_fill_ochre(palette="olsen_seq");SampUn
  scale_fill_manual(values = Trt)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        #axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_blank())+ 
  guides(fill=guide_legend(nrow = 2))+
  xlab("")+ylab("Relative Abundance");SampUn


Trt<-(SampUn + SampB + plot_layout(guides = 'keep'))+
  plot_annotation(title = 'Unburned                                            Burned',
   theme = theme(plot.title = element_text(size = 20),
           legend.position = "bottom"));Trt



#----
#----
#*********************************************************************************************************----
# --------------------------------- TREATMENT PER SAMPLE TYPE (2%)-------------------------------------------------------
#*********************************************************************************************************----


#At the pne percent level
RelGenBSt2$Genus[RelGenBSt2$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 
RelGenUSt2$Genus[RelGenUSt2$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 
RelGenUSt2$Genus[RelGenUSt2$Genus =="Burkholderia-Caballeronia-Paraburkholderia" ]<-"Burkholderia" 



Trt1<-c("< 2% abund."="#dee0df", "Acidibacter"="#c49993",
        "Bryobacter"="#baa4bf","Burkholderia"="#d4c9d6",
        "Can. Udaeobacter"="#faf7d2","Chitinophagaceae"="#c7c05f",
        "Chthoniobacter"="#827c2f","Clostridiales"="#e0d4b1",
        "Domibacillus"="#d6c07e","Gemmatimonas"="#d1d0b4",
        "Massilia"="#145405","Mucilaginibacter"="#63633b",
        "Mycobacterium"="#787612","Paenibacillus"="#425778",
        "RB41"="#96806f","Sphingomonas"="#d3e2e3",
        "Streptomyces"="#272829")


SampB1<-ggplot(data=RelGenBSt2, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack",show.legend = T) + #,show.legend = F
  #scale_fill_ochre(palette="olsen_seq");SampB1
  scale_fill_manual(values = Trt1)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        #axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
        legend.position="right",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=21))+ 
  guides(fill=guide_legend(ncol = 1))+
  xlab("")+ylab("");SampB1

SampUn1<-ggplot(data=RelGenUSt2, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack", show.legend = F) + 
  #scale_fill_ochre(palette="olsen_seq");SampUn1
  scale_fill_manual(values = Trt1)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        #axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
        legend.position="right",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=21))+ 
  guides(fill=guide_legend(ncol = 1))+
  xlab("")+ylab("Relative Abundance");SampUn1


Trt2<-(SampUn1 + SampB1 + plot_layout(guides = 'keep'))+
  plot_annotation(title = 'Unburned                           Burned',
         theme = theme(plot.title = element_text(size = 20),
          legend.position = "right"));Trt2




dir.create("Bacteria/Analysis/RelativeAbundance/Graphs/",recursive= T)
pdf("Bacteria/Analysis/RelativeAbundance/Graphs/Treatment-SampleType-Genus-Fun-3Per.pdf", height=6.5, width=12)
Trt
dev.off()


pdf("Bacteria/Analysis/RelativeAbundance/Graphs/Treatment-SampleType-Genus-Fun-2Per.pdf", height=6.5, width=12)
Trt2
dev.off()










#----
#----
#*********************************************************************************************************----
# --------------------------------- TIMES SINCE FIRES -------------------------------------------------------
#*********************************************************************************************************----
TsfU<-c("< 3% abund."="#dee0df","Anthelopsis"="#8a02a8","Balsamia"="#817399","Cenococcum"="#b2aabf",
        "Cladophialophora"="#534b61","Coniochaeta"="#6b2b13",
        "Cortinarius"="#351847","Geminibasidium"="#d1a82e", "Geopora"="#664c1c",
        "Hyaloscypha"="#453414","Hyphodiscus"="#c71636","Inocybe"="#B49A67",
        "Knufia"="#0207a8","Penicillium"="#607D8B","Pyronema" ="#324e5c",
        "Ramaria"="#0a0909","Sclerotinia"="#ebe2be","Thelephora"="#304f02",
        "Tomentella"="#193d11",
        "Trechispora"="#087527","Tubeufiales"="#9FB6C2",
        "Umbelopsis"="#9FB6C2","Venturia"="#0d1f00","unidentified"="#3d3d40")





TsfTU<-ggplot(data=RelGenTUtsf1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  ggtitle("Thawed")+
  #scale_fill_ochre(palette="olsen_seq")+
  scale_fill_manual(values = TsfU)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=18,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16,angle=90, colour = "black"),
        legend.position="right",legend.title=element_text(size=18),
        legend.text = element_text(face = "italic", size=16))+
  guides(fill=guide_legend(ncol = 1))+
  xlab("Time since fire (days)")+ylab("");TsfTU




TsfOU<-ggplot(data=RelGenOUtsf1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack", show.legend = F) + #,show.legend = F add after you see names and are rwady to plot
  #scale_fill_ochre(palette="olsen_seq")+ 
  ggtitle("Original")+
  scale_fill_manual(values = TsfU)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        legend.position="right",legend.title=element_text(size=18),
        legend.text = element_text(face = "italic", size=16))+
  guides(fill=guide_legend(ncol = 1))+ 
  xlab("Time since fire (days)")+ylab("");TsfOU


TsfFU<-ggplot(data=RelGenFUtsf1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack",show.legend = F) + 
  ggtitle("Frozen DNA")+
  #scale_fill_ochre(palette="olsen_seq")+
  scale_fill_manual(values = TsfU)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        legend.position="right",legend.title=element_text(size=18),
        legend.text = element_text(face = "italic", size=16))+ 
  xlab("Time since fire (days)")+ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));TsfFU
#guides(fill=guide_legend(nrow = 5));TsfFB



#----
#----
#*********************************************************************************************************----
# --------------------------------- BURNED SAMPLES -------------------------------------------------------
#*********************************************************************************************************----

TsfB<-c("< 3% abund."="#dee0df","Alternaria"="#946363","Amphinema"="#730030",
        "Aspergillus"="#3D0000", "Balsamia"="#817399",  "Bezerromyces"="#afbd19",
        "Coniochaeta"="#6b2b13","Coprinellus"="#8f4b07","Cortinarius"="#351847",
        "Fimetariella"="#876a35","Gelasinospora"="#5c4025","Geminibasidium"="#d1a82e",
        "Hyaloscypha"="#453414","Inocybe"="#B49A67",
        "Mycena"="#a86c6c",
       "Penicillium"="#607D8B",
        "Pyronema" ="#324e5c", "Tephrocybe" ="#4b6f42","Tomentella"="#193d11","Pyronemataceae"="#2b335c")
        
      


TsfTB<-ggplot(data=RelGenTBtsf1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack", show.legend = T) + 
  ggtitle("Thawed")+
  #scale_fill_ochre(palette="olsen_seq")+
  scale_fill_manual(values = TsfB)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=18,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        legend.position="right",legend.title=element_text(size=18),
        legend.text = element_text(face = "italic", size=16))+
        xlab("Time since fire (days)")+ylab("")+
        guides(fill=guide_legend(ncol = 1));TsfTB


TsfOB<-ggplot(data=RelGenOBtsf1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack",show.legend = F) + 
  ggtitle("Original")+
  #scale_fill_ochre(palette="olsen_seq")+
  scale_fill_manual(values = TsfB)+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.title=element_text(size=16),
        axis.title = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=16,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        legend.position="right",legend.title=element_text(size=18),
        legend.text = element_text(face = "italic", size=16))+ 
  xlab("Time since fire (days)")+ylab("")+
  guides(fill=guide_legend(ncol = 1));TsfOB


TsfFB<-ggplot(data=RelGenFBtsf1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack",show.legend = F) + 
  ggtitle("Frozen DNA")+
  #scale_fill_ochre(palette="olsen_seq")+
  scale_fill_manual(values = TsfB)+
  theme_bw() + theme(panel.grid = element_blank(),
      plot.title=element_text(size=16),
     axis.title = element_text(size = 18, colour = "black"), 
      axis.text.y = element_text(size=18,angle=90,hjust=0.60,colour = "black"),
      axis.text.x = element_text(size=16, colour = "black"),
      legend.position="right",legend.title=element_text(size=18),
     legend.text = element_text(face = "italic", size=16))+ 
  guides(fill=guide_legend(ncol =  1))+
  xlab("Time since fire (days)")+ylab("Relative Abundance");TsfFB
  #
#----
#*********************************************************************************************************----
# --------------------------------- CREATE PANEL OF PLOTS FOR EXPORT-----------------------------------------------
#*********************************************************************************************************----


BurnedPlots<-(TsfFB +TsfOB + TsfTB + plot_layout(guides = 'keep')) +
  plot_annotation(title = 'Burned',
   theme = theme(plot.title = element_text(size = 20)))+ 
   plot_layout(ncol = 3);BurnedPlots

UnburnPlots<-(TsfFU +TsfOU + TsfTU + plot_layout(guides = 'keep'))+
  plot_annotation(title = 'Unburned',
  theme = theme(plot.title = element_text(size = 20)))+ 
  plot_layout(ncol = 3);UnburnPlots

Trt<-(SampUn + SampB + plot_layout(guides = 'keep'))+
  plot_annotation(title = '',
  theme = theme(plot.title = element_text(size = 20)))+ 
  plot_layout(nrow = 1);Trt


#Export individual plots before making panels ...........................................
dir.create(("Analysis/RelativeAbundance/Graphs"), recursive=TRUE)

pdf("Analysis/RelativeAbundance/Graphs/Burned-TSF-SampleType-Genus-Fun.pdf", height=8, width=12)
BurnedPlots
dev.off()

pdf("Analysis/RelativeAbundance/Graphs/Unburned-TSF-SampleType-Genus-Fun.pdf", height=8, width=12)
UnburnPlots
dev.off()

pdf("Analysis/RelativeAbundance/Graphs/Treatment-SampleType-Genus-Fun.pdf", height=6.5, width=12)
Trt
dev.off()



pdf("Analysis/RelativeAbundance/ThawedSample-Contam.pdf", height=8, width=9)
TsfTU
dev.off()

