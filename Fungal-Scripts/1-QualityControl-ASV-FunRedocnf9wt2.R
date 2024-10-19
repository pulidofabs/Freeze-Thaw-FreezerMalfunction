#Reset R's Brain
rm(list=ls())

#Set working directory-----------------------------------------------------------------------------------
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/NoCNF9WT2/")

#Load new packages -------------------------------------------------------
library(tidyverse) #required to load tax table
library(dplyr)
library(EcolUtils) 
library(SPECIES)
library(BiodiversityR)
library(scales)
library(plyr)


#LOAD DATA ------------------------------------------------------------------------------------------------------
Metadata<-read.csv("Fungi/Metadata/Fungi-MergedAll-NC.csv", row.names = 1,header = TRUE)#make sure to load metadata w/o negative controls
dim(Metadata)#57x13

# * * LOAD OTU TABLES (no singletons, contm'd samples removed-------------------------------------------------------------------------------
RawOtu<-read.csv("Fungi/ASVtables/Fungi-Table-FilterTaxon-1-4RT-9WT2.csv",row.names = 1, check.names = FALSE)
dim(RawOtu)#3193x60


#*********************************************************************************************************************************------
#-------------CREATE CLEAN ASV TABLE AND ID/TAXA LABEL TABLES  ---------------------------------------------
#*********************************************************************************************************************************------
# * Remove FeatureID and Taxonomy Labels to use later................................................................................
lastcolumn <- ncol(RawOtu); lastcolumn #--60---What is the last column in the dataset, need for removing taxa
taxononomy<-RawOtu[,60]; taxononomy #Extract taxonomy column from the dataset
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
RawOtuTable  <- RawOtu2[ ,1:59];head(RawOtuTable)#number as shown above-1

# TRANSPOSE RAW OTU TABLE............................................................................................
OtuTrans <- t(RawOtuTable);dim(OtuTrans)#59x3193



#EXPORT DATA -------------------------------------------------------------------------------------------------------
dir.create(file.path("Fungi/ASVtables"), recursive = TRUE)

write.csv(RawOtuTable, "Fungi/ASVtables/RawOtuTable.csv")
write.csv(OtuTrans, "Fungi/ASVtables/OtuTrans.csv")
write.csv(TaxonID, "Fungi/ASVtables/TaxononomyLabels.csv")


#----
#----
#*********************************************************************************************************************************------
# --- -----------   QUALITY CONTROL  ---------------------------------------------
#*********************************************************************************************************************************------
colnames(RawOtuTable)#Look at rownames to see the name of neg + pos controls

# * look at Negative & Mocl DNA controls ................................................................
#NegControls <- which(colnames(RawOtuTable) %in% c("")); NegControls# manually removed
MockComm<-which(colnames(RawOtuTable) %in%c("PositivePlate1_F","FPCRPositiveCtrl"));MockComm#27x59

# make otu table of Negative DNA controls & Mock Community...............................................
OtuMockComm<-RawOtuTable[,MockComm]; head(OtuMockComm)
RawOtuTable<-RawOtuTable[,-MockComm]; dim(RawOtuTable)#3193x57

#Remove Negcontrols and Moc comunities from the RawOtu Table.............................................
#OtuNegControl<-RawOtuTable[ ,NegControls]; head(OtuNegControl)
#RawOtuTable<-RawOtuTable[,-NegControls]; dim(RawOtuTable)#5039X13

colnames(RawOtuTable)#verify removal, if not remove, rerun the Mock scripts & it removes

#export table w out controls ...........................................................................
write.csv(RawOtuTable, "Fungi/ASVtables/RawOtuTableNC-CNF9WT2.csv")

#Transpose table to put it in correct order for analysis................................................
OtuTransNC <- t(RawOtuTable);dim(OtuTransNC) #57x3193
write.csv(OtuTransNC, "Fungi/ASVtables/OtuTransNC-CNF9WT2.csv")


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
OtuTrans2<-OtuTransNC[rowSums(OtuTransNC)>5786,] 
dim(OtuTrans2)#57x3193

# * * Remove low abundance seq, normalize to 13911 seq/sample, no iterations (1x)...................
OtuRare1 <- rrarefy(OtuTrans2, 5786)
dim(OtuRare1)#57x57x3193

# *  * Normalize the data 100x and get the mean.....................................................
OtuRare2<- rrarefy.perm(OtuTrans2, sample = 5786, n = 250, round.out = T)
dim(OtuRare2)#57x3193


# * * * TEST IF RAREFACTION METHOD MATTERS..........................................................
mantel(OtuRare1, OtuRare2)#0.9858

# REMOVE COLUMNS WITH ZERO READS (OTURARE2)........................................................
# * Create an object containing zeros .........................
zeroes <- which(colSums(OtuRare2)==0)
head(zeroes)


# * Remove the columns containing zero's from OTU table ........
OtuRare3NC <- OtuRare2[ , -zeroes]
head(OtuRare3NC[,1:2]);dim(OtuRare3NC)#57X2647


# * * EXPORT RAREFACTION TABLES  ......................................................
write.csv(OtuTrans2, "Fungi/ASVtables/OtuTrans2-CNF9WT2.csv")
write.csv(OtuRare2, "Fungi/ASVtables/OtuRare2-CNF9WT2.csv")
write.csv(OtuRare1, "Fungi/ASVtables/OtuRare1-CNF9WT2.csv")
write.csv(OtuRare3NC, "Fungi/ASVtables/OtuRare3NZ-CNF9WT2.csv")#w/o zeroes



#----
#----
#***************************************************************************************************************************************----
#RAREFY TABLE TO MATCH THE RAREFIED TABLE L-------------------------------------------------------------------------------------------------
#***************************************************************************************************************************************----
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3
metadata<-read.csv("Fungi/Metadata/Fungi-MergedAll-NC.csv", na.strings = "N/A", header = TRUE);dim(metadata)#57x14
OtuRare3NC<-read.csv("Fungi/ASVtables/OtuRare3NZ-CNF9WT2.csv", check.names = FALSE);dim(OtuRare3NC)#57x2648

#Verify rownames match 
row.names(OtuRare3NC)==row.names(metadata)#yes and reimport w/o rownames


names(OtuRare3NC)[1]<- "SampleID"
names(metadata)[1]<- "SampleID"

#Verify that column1 has the same name in both tables..............................................
names(metadata[,1:2]); names(OtuRare3NC[,1:2])
dim(metadata);dim(OtuRare3NC)


#Rarefy Metadata to match rarefied ASV table........................................................
MetaRare<-metadata %>% semi_join(OtuRare3NC, by = "SampleID") # keep rows with matching ID
dim(MetaRare)#57x14


#Export rarefied metadata .........................................................................
dir.create("Fungi/Metadata")
write.csv(MetaRare, "Fungi/Metadata/MetaRareFun-CNF9WT2.csv")



#----
#----
#***************************************************************************************************************************************----
# CALCULATE ALPHA DIVERSITY -------------------------------------------------------------------------------------------------
#***************************************************************************************************************************************----
# RICHNESS METHOD 1: BIODIVERSITY CALCULATED-........................................................................

# * * *  Load metadata that correspands to R-rarefied table created above .....
MetaRare<-read.csv("Fungi/Metadata/MetaRareFun-CNF9WT2.csv", row.names = 1);dim(MetaRare)#57x13
OtuRare3NC<-read.csv("Fungi/ASVtables/OtuRare3NZ-CNF9WT2.csv", row.names = 1,check.names = FALSE);dim(OtuRare3NC)#57x2647

#Verify that rownames match
all(rownames(MetaRare)==rownames(OtuRare3NC))#Yes


# * Richness calculations using Biodivesity package-----------------------------------
OtuRichness <- estimateR(OtuRare3NC)
OtuRichness<- t(estimateR(OtuRare3NC))#run function on transformed table to switch rows and columns

# CREATE FUNCTION TO ESTIMATE DIVERSITY METRICS ------------------------------
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
dim(OtuSppRichness)#57X2655

# * * * Export Biodiversity package spp richness--------------------------------------------------------------------
dir.create(file.path("Fungi/Analysis/Diversity/Tables"), recursive = TRUE)
write.csv(OtuSppRichness, "Fungi/Analysis/Diversity/Tables/CoreMetricsAll-CNF9WT2.csv")

#sUBSET ONLY THE RICHNESS METRICS
ncol(OtuSppRichness)#2655

#1 column and 6th column
CoreMetrics1<-OtuSppRichness[,2649:2655];head(CoreMetrics1)#2655-6
CoreMetrics2<-as.data.frame(OtuSppRichness[,1]);head(CoreMetrics2)
names(CoreMetrics2)[1]<- "ShannonH"

CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)

#Export core metrics results.............................................................
write.csv(CoreMetrics,"Fungi/Analysis/Diversity/Tables/AlphaCoreMetrics-CNF9WT2.csv")

#----
#----
#****************************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#****************************************************************************************************************----
#Importa tables with row.names in row.names=1 jsut to verify if names match and then reimport without row.names
##for analysis

#Merged table manually code seems to be off
CoreMetrics<-read.csv("Fungi/Analysis/Diversity/Tables/AlphaCoreMetrics-CNF9WT2.csv", check.names = FALSE)
MetaRare<-read.csv("Fungi/Metadata/MetaRareFun-CNF9WT2.csv");dim(MetaRare)#57x14
names(CoreMetrics)[1]<- "SampleID"

all(row.names(CoreMetrics)==row.names(MetaRare))

#Add the coremetrics file to the metadata file................................
MetaRareSpp<-merge(MetaRare,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#57x22

#Export rarefied metadata w core metrics attached-------------------------------
write.csv(MetaRareSpp, "Fungi/Metadata/MetaRareSppFun-CNF9WT2.csv")

#..................................................................................
MetaRareSpp<-read.csv("Fungi/Metadata/MetaRareSppFun-CNF9WT2.csv")

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

write.csv(MetaSppOrig, "Fungi/Metadata/MetaSppOrig-CNF9WT2.csv")
write.csv(MetaSppThaw, "Fungi/Metadata/MetaSppRaw-CNF9WT2.csv")
write.csv(MetaSppFrozen, "Fungi/Metadata/MetaSppFrozen-CNF9WT2.csv")
