#Reset R's Brain
rm(list=ls())

#Set working directory-----------------------------------------------------------------------------------
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/Lib7Thawed-NoCNF9WT2")

#Load new packages -------------------------------------------------------
library(tidyverse) #required to load tax table
library(dplyr)
library(EcolUtils) 
library(SPECIES)
library(BiodiversityR)
library(scales)
library(plyr)


#LOAD DATA ------------------------------------------------------------------------------------------------------
Metadata<-read.csv("Bacteria/Metadata/Bacteria-MergedAll.csv", row.names = 1,header = TRUE)#make sure to load metadata w/o negative controls
dim(Metadata)#56x14

# * * LOAD OTU TABLES (no singletons, contm'd samples removed-------------------------------------------------------------------------------
RawOtu<-read.csv("Qiime/Merged/Bacteria/ExportedRaw/Bacteria-Table-FilterTaxon-1-4RT-9WT2.csv",row.names = 1, check.names = FALSE)
dim(RawOtu)#8691x626


#*********************************************************************************************************************************------
#-------------CREATE CLEAN ASV TABLE AND ID/TAXA LABEL TABLES  ---------------------------------------------
#*********************************************************************************************************************************------
# * Remove FeatureID and Taxonomy Labels to use later................................................................................
lastcolumn <- ncol(RawOtu); lastcolumn #--62---What is the last column in the dataset, need for removing taxa
taxononomy<-RawOtu[,62]; taxononomy #Extract taxonomy column from the dataset
featureID <- row.names(RawOtu); featureID #Extract feature ID, from the dataset
TaxonID <- data.frame(featureID,taxononomy); head(TaxonID) #Create a dataframe of the featureID and taxonomy

# * Clean Taxonomy labels by separating them by (;)..................................................................................
   #divide last row into subsets #separate taxonomy w spaces 
SplitTaxonomy<-ldply(str_split(string = TaxonID$taxononomy, pattern=";"), rbind)#Separate columns w ";" & convert to dataframe
names(SplitTaxonomy)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #rename the column names
SplitTaxonomy2 <- as.data.frame(lapply(SplitTaxonomy, gsub, pattern=" ", replacement=""));SplitTaxonomy2
tail(SplitTaxonomy2);head(SplitTaxonomy2)

# * Change the format of the TaxonID dataframe created above so that is is split.....................................................
head(TaxonID)
TaxonID <- cbind(TaxonID[,1:2 ], SplitTaxonomy2);head(TaxonID)
rownames(TaxonID) <- featureID;head(TaxonID)#add feaure ID to table

# * * Look at your dataset to make sure that they have the same length...............................................
dim(TaxonID);dim(RawOtu) # should have same row length 

# * ATTACH TAXON ID CHANGES TO RAWOTU TABLE TO CREATE A CLEAN TABLE
RawOtu2<- cbind(RawOtu, TaxonID);head(RawOtu2[1:2,])

# CREATE SV TABLE--NO TAXONOMY ......................................................................................
RawOtuTable  <- RawOtu2[ ,1:61];head(RawOtuTable)#number as shown above-1

# TRANSPOSE RAW OTU TABLE............................................................................................
OtuTrans <- t(RawOtuTable);dim(OtuTrans)#61x8691



#EXPORT DATA -------------------------------------------------------------------------------------------------------
dir.create(file.path("Bacteria/ASVtables"), recursive = TRUE)

write.csv(RawOtuTable, "Bacteria/ASVtables/RawOtuTable.csv")
write.csv(OtuTrans, "Bacteria/ASVtables/OtuTrans.csv")
write.csv(TaxonID, "Bacteria/ASVtables/TaxononomyLabels.csv")


#----
#----
#*********************************************************************************************************************************------
# --- -----------   QUALITY CONTROL  ---------------------------------------------
#*********************************************************************************************************************************------
colnames(RawOtuTable)#Look at rownames to see the name of neg + pos controls

# * look at Negative & Mocl DNA controls ................................................................
#Make out table of controsl and then remove them..........
MockComm<-which(colnames(RawOtuTable) %in%c("PositivePlate1_B","BPCRPositiveCtrl"));MockComm#27x59
OtuMockComm<-RawOtuTable[,MockComm]; dim(RawOtuTable)#8691x61
RawOtuTable<-RawOtuTable[,-MockComm]; dim(RawOtuTable)#8691x59

#Remove Negcontrols and Moc comunities from the RawOtu Table.............................................
NegControls <- which(colnames(RawOtuTable) %in% c("negative_Plate2_B","NegativePlate1_B","BPCRNegCtrl"));NegControls# manually removed
OtuNegComm<-RawOtuTable[,NegControls]; dim(RawOtuTable)#8691x59
RawOtuTable<-RawOtuTable[,-NegControls]; dim(RawOtuTable)#8691x56

colnames(RawOtuTable)#verify removal, if not remove, rerun the Mock scripts & it removes

#export table w out controls ...........................................................................
write.csv(RawOtuTable, "Bacteria/ASVtables/RawOtuTableNC-CNF9WT2.csv")

#Transpose table to put it in correct order for analysis................................................
OtuTransNC <- t(RawOtuTable);dim(OtuTransNC) #56x8691
write.csv(OtuTransNC, "Bacteria/ASVtables/OtuTransNC-CNF9WT2.csv")


#verify that row names match with metadata .............................................................
all(row.names(OtuTransNC)==row.names(Metadata))


#Look at total sequences/sample to select rarefactin depth...............................................
TotSeq<-sort(rowSums(OtuTransNC), decreasing = TRUE);TotSeq #lowest is 5787 

dir.create(file.path("Analysis/QualityControl"), recursive = T)
write.csv(TotSeq, "Analysis/QualityControl/TotalSequence-Samples-CNF9WT2.csv")



#.----
#.----
#***********************************************************************************************************************------------
# RAREFY TABLE- FULL Bacteria--------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************------------
# Chose rarefaction value using qiime feature table, removed botto 10%

# * * Trial 1 Remove samples w not enough sequences.....................................................................
# * * * Subset data to maintain only rows with totals above rarefaction depth
OtuTrans2<-OtuTransNC[rowSums(OtuTransNC)>6493,] 
dim(OtuTrans2)#56x8691

# * * Remove low abundance seq, normalize to 13911 seq/sample, no iterations (1x)...................
OtuRare1 <- rrarefy(OtuTrans2, 6493)
dim(OtuRare1)#56x8691

# *  * Normalize the data 100x and get the mean.....................................................
OtuRare2<- rrarefy.perm(OtuTrans2, sample = 6493, n = 250, round.out = T)
dim(OtuRare2)#56x8691


# * * * TEST IF RAREFACTION METHOD MATTERS..........................................................
mantel(OtuRare1, OtuRare2)#0.9858

# REMOVE COLUMNS WITH ZERO READS (OTURARE2)........................................................
# * Create an object containing zeros .........................
zeroes <- which(colSums(OtuRare2)==0)
head(zeroes)


# * Remove the columns containing zero's from OTU table ........
OtuRare3NC <- OtuRare2[ , -zeroes]
head(OtuRare3NC[,1:2]);dim(OtuRare3NC)#56x8086


# * * EXPORT RAREFACTION TABLES  ......................................................
write.csv(OtuTrans2, "Bacteria/ASVtables/OtuTrans2-CNF9WT2.csv")
write.csv(OtuRare2, "Bacteria/ASVtables/OtuRare2-CNF9WT2.csv")
write.csv(OtuRare1, "Bacteria/ASVtables/OtuRare1-CNF9WT2.csv")
write.csv(OtuRare3NC, "Bacteria/ASVtables/OtuRare3NZ-CNF9WT2.csv")#w/o zeroes



#----
#----
#***************************************************************************************************************************************----
#RAREFY TABLE TO MATCH THE RAREFIED TABLE L-------------------------------------------------------------------------------------------------
#***************************************************************************************************************************************----
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3
metadata<-read.csv("Bacteria/Metadata/Bacteria-MergedAll.csv", na.strings = "N/A", header = TRUE);dim(metadata)#56x15
OtuRare3NC<-read.csv("Bacteria/ASVtables/OtuRare3NZ-CNF9WT2.csv", check.names = FALSE);dim(OtuRare3NC)#56x8087

#Verify rownames match 
row.names(OtuRare3NC)==row.names(metadata)#yes and reimport w/o rownames

names(OtuRare3NC)[1]<- "SampleID"
names(metadata)[1]<- "SampleID"

#Verify that column1 has the same name in both tables..............................................
names(metadata[,1:2]); names(OtuRare3NC[,1:2])
dim(metadata);dim(OtuRare3NC)


#Rarefy Metadata to match rarefied ASV table........................................................
MetaRare<-metadata %>% semi_join(OtuRare3NC, by = "SampleID") # keep rows with matching ID
dim(MetaRare)#56x15


#Export rarefied metadata .........................................................................
dir.create("Bacteria/Metadata")
write.csv(MetaRare, "Bacteria/Metadata/MetaRareBac-CNF9WT2.csv")



#----
#----
#***************************************************************************************************************************************----
# CALCULATE ALPHA DIVERSITY -------------------------------------------------------------------------------------------------
#***************************************************************************************************************************************----
# RICHNESS METHOD 1: BIODIVERSITY CALCULATED-........................................................................

# * * *  Load metadata that correspands to R-rarefied table created above .....
MetaRare<-read.csv("Bacteria/Metadata/MetaRareBac-CNF9WT2.csv", row.names = 1);dim(MetaRare)#56x15
OtuRare3NC<-read.csv("Bacteria/ASVtables/OtuRare3NZ-CNF9WT2.csv", row.names = 1,check.names = FALSE);dim(OtuRare3NC)#56x8083

#Verify that rownames match
all(rownames(MetaRare)==rownames(OtuRare3NC))#Yes


# * Richness calculations using Biodivesity package-----------------------------------
OtuRichness <- estimateR(OtuRare3NC)
OtuRichness<- t(estimateR(OtuRare3NC))#run Bacction on transformed table to switch rows and columns

# CREATE BacCTION TO ESTIMATE DIVERSITY METRICS ------------------------------
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  # pdf(paste("figures/richnesscores_",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", 
       ylab="chao", col=alpha("red", 0.5),pch=16)
  mtext("A",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",
       ylab="ACE",col=alpha("black", 0.5),pch=16)
  mtext("B",side=3,adj=0)
  #dev.off()
}

# DIVERSITY INDICES ------------------------------------------------------------------------------------------------
shanEntro <- diversity(OtuRare3NC) # Shannon entropy
shannon <- exp(OtuRare3NC) ## Shannon number of diversity
simpson <- diversity(OtuRare3NC, "simpson")#Simpson diversity
simpEven<- diversity(OtuRare3NC, "inv") ## Simpson Inverse (Eveness)

# * Dataframe of shannon entropy, diversity &  simpson diversity--------------
otu.richness <- data.frame(shanEntro, shannon, simpson, simpEven)

# * Add above diversity metrics to S obs, chao1, ACE
OtuSppRichness <- cbind(otu.richness, OtuRichness)
dim(OtuSppRichness)#56X8091

# * * * Export Biodiversity package spp richness--------------------------------------------------------------------
dir.create(file.path("Bacteria/Analysis/Diversity/Tables"), recursive = TRUE)
write.csv(OtuSppRichness, "Bacteria/Analysis/Diversity/Tables/CoreMetricsAll-CNF9WT2.csv")

#sUBSET ONLY THE RICHNESS METRICS
ncol(OtuSppRichness)#8091

#1 column and 6th column
CoreMetrics1<-OtuSppRichness[,8085:8091];head(CoreMetrics1)#8091-6
CoreMetrics2<-as.data.frame(OtuSppRichness[,1]);head(CoreMetrics2)
names(CoreMetrics2)[1]<- "ShannonH"

CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)

#Export core metrics results.............................................................
write.csv(CoreMetrics,"Bacteria/Analysis/Diversity/Tables/AlphaCoreMetrics-CNF9WT2.csv")

#----
#----
#****************************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#****************************************************************************************************************----
#Importa tables with row.names in row.names=1 jsut to verify if names match and then reimport without row.names
##for analysis

#Merged table manually code seems to be off
CoreMetrics<-read.csv("Bacteria/Analysis/Diversity/Tables/AlphaCoreMetrics-CNF9WT2.csv", check.names = FALSE)
MetaRare<-read.csv("Bacteria/Metadata/MetaRareBac-CNF9WT2.csv");dim(MetaRare)#56x15
names(CoreMetrics)[1]<- "SampleID"

all(row.names(CoreMetrics)==row.names(MetaRare))#True

#Add the coremetrics file to the metadata file................................
MetaRareSpp<-merge(MetaRare,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#56x23

#Export rarefied metadata w core metrics attached-------------------------------
write.csv(MetaRareSpp, "Bacteria/Metadata/MetaRareSppBac-CNF9WT2.csv")

#..................................................................................
MetaRareSpp<-read.csv("Bacteria/Metadata/MetaRareSppBac-CNF9WT2.csv")

#Separate the metadata --Thawed vs Original.........................................
attach(MetaRareSpp)
MetaSppOrig<-MetaRareSpp[which(SampleType== "Original"), ];head(MetaSppOrig[,1:9]);dim(MetaSppOrig)#19x23
MetaSppThaw<-MetaRareSpp[which(SampleType == "Thawed"), ];head(MetaSppThaw[,1:9]);dim(MetaSppThaw)#19x23
MetaSppFrozen<-MetaRareSpp[which(SampleType == "FrozenDNA"), ];head(MetaSppFrozen[,1:9]);dim(MetaSppFrozen)#19x23
detach(MetaRareSpp)

#Combine tables using SampleName.....................................................
#MetaRareSppCor<-MetaSppOrig %>% join(MetaSppThaw, by = "SampleName") # keep rows with matching ID
#dim(MetaRareSppCor);head(MetaRareSppCor[,1:30])#14x9

#Export table and remove unecessary column...........................................................
#write.csv(MetaRareSppCor, "Analysis/Metadata/MetaRareSppCor.csv")

write.csv(MetaSppOrig, "Bacteria/Metadata/MetaSppOrig-CNF9WT2.csv")
write.csv(MetaSppThaw, "Bacteria/Metadata/MetaSppRaw-CNF9WT2.csv")
write.csv(MetaSppFrozen, "Bacteria/Metadata/MetaSppFrozen-CNF9WT2.csv")

