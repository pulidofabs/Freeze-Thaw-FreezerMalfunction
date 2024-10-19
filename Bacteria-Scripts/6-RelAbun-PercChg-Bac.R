#Reset R's Brain
rm(list=ls())


#Set working directory................................................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/Lib7Thawed-NoCNF9WT2/")



#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)



#Load data--------------------------------------------------------------------------------------------------------
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


#--Subset data by treatment (burned vs unburned) ......................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) 

# * * Subset All to maintain only the burned samples of all data......................
physeqB<-subset_samples(physeq,Treatment=="Burned");physeqB#7926x30
physeqUn<-subset_samples(physeq,Treatment=="Unburned");physeqUn#7926x24

#Set base level to original ............................................................




#.----
#********************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -----------------------------------------------------------------
#*********************************************************************************************-----

#* * Burned ........................................................................
GenBB<- tax_glom(physeqB, "Genus")
RelGenBB <- transform_sample_counts(GenBB, function(x) x / sum(x))
RelGenBB1<-psmelt(RelGenBB)

#* * Unburned ......................................................................
GenBU<- tax_glom(physeqUn, "Genus")
RelGenBU<- transform_sample_counts(GenBU, function(x) x / sum(x))
RelGenBU1<-psmelt(RelGenBU)


#Quality control ..................................................................
RelGenBB1$SampleType<-as.factor(RelGenBB1$SampleType)
RelGenBU1$SampleType<-as.factor(RelGenBU1$SampleType)

RelGenBB1$SampleType<-try(relevel(RelGenBB1$SampleType,ref="Original"))
RelGenBU1$SampleType<-try(relevel(RelGenBU1$SampleType,ref="Original"))


#.----
#.----
#**************************************************************************************************************************-----
# EXPORT DATA FILES --------------------------------------------------------------------------------------------------------------
#**************************************************************************************************************************-----
dir.create(file.path("Bacteria/Analysis/RelativeAbundance/Tables/DesStats"), recursive=TRUE)


#Re-Import soil samples ......................................................................................................
write.csv(RelGenBB1,"Bacteria/Analysis/RelativeAbundance/Tables/DesStats/RelAbund-Burn-All.csv")
write.csv(RelGenBU1,"Bacteria/Analysis/RelativeAbundance/Tables/DesStats/RelAbund-Unburn-All.csv")


#---
#-----
#********************************************************************************************************************----
#---DESCRIPTIVE STATISTICS -------------------------------------------------------------------------------------------------
#********************************************************************************************************************----
RelGenBB1<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/DesStats/RelAbund-Burn-All.csv")
RelGenBU1<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/DesStats/RelAbund-Unburn-All.csv")


#Subset data only thawed and original 
BacBOT <- subset(RelGenBB1, SampleType == "Original" | SampleType == "Thawed");unique(BacBOT$SampleType)
BacBOF <- subset(RelGenBB1, SampleType == "Original" | SampleType == "Frozen");unique(BacBOF$SampleType)

BacUOT <- subset(RelGenBU1, SampleType == "Original" | SampleType == "Thawed");unique(BacUOT$SampleType)
BacUOF <- subset(RelGenBU1, SampleType == "Original" | SampleType == "Frozen");unique(BacUOF$SampleType)


#Relevel to make original the base level
BacBOT$SampleType<-try(relevel(BacBOT$SampleType,ref="Original"))
BacBOF$SampleType<-try(relevel(BacBOF$SampleType,ref="Original"))

BacUOT$SampleType<-try(relevel(BacUOT$SampleType,ref="Original"))
BacUOF$SampleType<-try(relevel(BacUOF$SampleType,ref="Original"))

#Calculate percent changes....................................................................................----
library(dplyr)

#Burned orig-thaw
BBOT<-BacBOT %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();head(BBOT)

#Orig-Frozen
BBOF<-BacBOF %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();head(BBOF)


#dir_path <- "Bacteria/Analysis/RelativeAbundance/Tables/DeStats/"
#if (!file.exists(dir_path)) {
#  dir.create(dir_path, recursive = TRUE)
#}
write.csv(BBOT,"Bacteria/Analysis/RelativeAbundance/Tables/DesStats/Bac-B-OrigThaw-PerChange-All.csv")
write.csv(BBOF,"Bacteria/Analysis/RelativeAbundance/Tables/DesStats/Bac-B-OrigFrozen-PerChange-All.csv")


#Same stats but for Unburned dataset
#Burned orig-thaw
BUOT<-BacUOT %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();BUOT

#Orig-Frozen
BUOF<-BacUOF %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();BUOF



write.csv(BUOT,"Bacteria/Analysis/RelativeAbundance/Tables/DesStats/Bac-U-OrigThaw-PerChange-All.csv")
write.csv(BUOF,"Bacteria/Analysis/RelativeAbundance/Tables/DesStats/Bac-U-OrigFrozen-PerChange-All.csv")






#----
#-----
#***********************************************************************************************************************-----
# CALCULATE SIGNIFICANCE----------------------------------------------------------------------------------------------------
#***********************************************************************************************************************-----
#Re-Import soil samples ......................................................................................................
RelGenBB1<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/DesStats/RelAbund-Burn-All.csv")
RelGenBU1<-read.csv("Bacteria/Analysis/RelativeAbundance/Tables/DesStats/RelAbund-Unburn-All.csv")



#Burned-----------------------------------------------------------------------------------------------
attach(RelGenBB1)
GemmatimonasB<-RelGenBB1[which(RelGenBB1$Genus == "Gemmatimonas" ),]; head(GemmatimonasB$Genus)
DomibacillusB<-RelGenBB1[which(RelGenBB1$Genus == "Domibacillus" ),]; head(DomibacillusB$Genus)
MassiliaB<-RelGenBB1[which(RelGenBB1$Genus == "Massilia" ),]; head(MassiliaB$Genus)
PaenibacillusB<-RelGenBB1[which(RelGenBB1$Genus == "Paenibacillus"),]; head(PaenibacillusB$Genus)
RB41B<-RelGenBB1[which(RelGenBB1$Genus == "RB41"),]; head(RB41B$Genus)
ClostridialesB<-RelGenBB1[which(RelGenBB1$Genus == "Clostridiales-Family XVIII" ),]; head(ClostridialesB$Genus)


#CHECK DATA FOR NORMALITY..........................................................................
"Burkholderia-Caballeronia-Paraburkholderia"

#Check histogram and run shapiro test ...........................................
hist(Bryobacter$Abundance);shapiro.test(Bryobacter$Abundance)
hist(Mycobacterium$Abundance);shapiro.test(Mycobacterium$Abundance)#normal
hist(Burkholderia$Abundance);shapiro.test(Burkholderia$Abundance)#normal



#-----
#----
#*****************************************************************************************************************************************----
#SIGNIFICANCE PER TAXA PER TREATMENT----------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************----
library(MASS)
library(lmerTest)
library(lme4)
library(MuMIn)
library(car)#anova for glmer
library(multcomp) #for post-hoc multilevel test#4


#Burned samples 
ASVneg1 <- lmer(Abundance ~ SampleType + (1|Plot)+(1|TSFdays)+(1|TSFdays), data = MassiliaB)#model failed to converge
ASVneg2 <- lmer(Abundance ~ SampleType  +  (1|Plot)+(1|TSFdays), data = MassiliaB)
ASVneg3 <- lmer(Abundance ~ SampleType  + (1|Subplot)+(1|TSFdays), data = MassiliaB)
ASVneg4 <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = MassiliaB)
ASVneg5 <- lmer(Abundance ~ SampleType  +  (1|Plot), data = MassiliaB)
ASVneg6 <- lmer(Abundance ~ SampleType  + (1|Subplot), data = MassiliaB)

#Check model to see which is better
AICc(ASVneg1,ASVneg2,ASVneg3,ASVneg4,ASVneg5,ASVneg6)#ASVneg5 better



#Quality control, set base level to original
#Burned samples ......................................................................
GemmatimonasB$SampleType<-as.factor(GemmatimonasB$SampleType)
DomibacillusB$SampleType<-as.factor(DomibacillusB$SampleType)
MassiliaB$SampleType<-as.factor(MassiliaB$SampleType)
PaenibacillusB$SampleType<-as.factor(PaenibacillusB$SampleType)
RB41B$SampleType<-as.factor(RB41B$SampleType)
ClostridialesB$SampleType<-as.factor(ClostridialesB$SampleType)

GemmatimonasB$SampleType<-try(relevel(GemmatimonasB$SampleType,ref="Original"))
DomibacillusB$SampleType<-try(relevel(DomibacillusB$SampleType,ref="Original"))
MassiliaB$SampleType<-try(relevel(MassiliaB$SampleType,ref="Original"))
PaenibacillusB$SampleType<-try(relevel(PaenibacillusB$SampleType,ref="Original"))
RB41B$SampleType<-try(relevel(RB41B$SampleType,ref="Original"))
ClostridialesB$SampleType<-try(relevel(ClostridialesB$SampleType,ref="Original"))



#Test for significance-..............................................................................................
#Massilia burned best is plot only
GemmaLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = GemmatimonasB)
Anova(GemmaLm, type=3)#0.2646
#PHGemma<-summary(glht(GemmaLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHGemma#T-O (0.000499) T-F(0.005632)

DomiLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = DomibacillusB)
Anova(DomiLm, type=3)#0.22400  
#PHDomi<-summary(glht(DomiLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHDomi#T-O (0.000499) T-F(0.005632)

MassLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = MassiliaB)
Anova(MassLm, type=3)#0.802223 
#PHMass<-summary(glht(MassLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHMass#T-O (0.000499) T-F(0.005632)

PaeniLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = PaenibacillusB)
Anova(PaeniLm, type=3)# 0.569613 
#PHPaeni<-summary(glht(PaeniLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHPaeni#T-O (0.000499) T-F(0.005632)

RB41BLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = RB41B)
Anova(RB41BLm, type=3)#0.96994 
#PHRB41B<-summary(glht(RB41BLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHRB41B#T-O (0.000499) T-F(0.005632)


ClostLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = ClostridialesB)
Anova(ClostLm, type=3)#0.0610 
#PHClost<-summary(glht(ClostLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHClost#T-O (0.000499) T-F(0.005632)




# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt"), recursive=TRUE)

#Export ANOVA values ...............................................................................................................
capture.output(Anova(Bryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Bryobacterglmer.csv")
capture.output(Anova(Bryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Can.Udaeobacterglmer.csv")
capture.output(Anova(Udaeobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Sphingomonasglmer.csv")


#-----
#-------
#************************************************************************************************************************----
# -- UNBURNED COMMUNITIES----------------------------------------------------------
##***********************************************************************************************************************----

#Unburned data--------------------------------------------------------------------------------------------------------
Acidibacter<-RelGenBU1[which(RelGenBU1$Genus == "Acidibacter"),]; head(Acidibacter$Genus)
Bryobacter<-RelGenBU1[which(RelGenBU1$Genus == "Bryobacter"),]; head(Bryobacter$Genus)
Burkholderia<-RelGenBU1[which(RelGenBU1$Genus == "Burkholderia-Caballeronia-Paraburkholderia"),]; head(Burkholderia$Genus)
CanUndaeobacter<-RelGenBU1[which(RelGenBU1$Genus == "Candidatus Udaeobacter" ),]; head(CanUndaeobacter$Genus)
Chitinophagaceae<-RelGenBU1[which(RelGenBU1$Genus == "Chitinophagaceae" ),]; head(Chitinophagaceae$Genus)
Chthoniobacter<-RelGenBU1[which(RelGenBU1$Genus == "Chthoniobacter" ),]; head(Chthoniobacter$Genus)
Gemmatimonas<-RelGenBU1[which(RelGenBU1$Genus == "Gemmatimonas" ),]; head(Gemmatimonas$Genus)
Massilia<-RelGenBU1[which(RelGenBU1$Genus == "Massilia" ),]; head(Massilia$Genus)
Mucilaginibacter<-RelGenBU1[which(RelGenBU1$Genus == "Mucilaginibacter" ),]; head(Mucilaginibacter$Genus)
Mycobacterium<-RelGenBU1[which(RelGenBU1$Genus == "Mycobacterium" ),]; head(Mycobacterium$Genus)
Paenibacillus<-RelGenBU1[which(RelGenBU1$Genus == "Paenibacillus" ),]; head(Paenibacillus$Genus)
RB41<-RelGenBU1[which(RelGenBU1$Genus == "RB41"),]; head(RB41$Genus)
Sphingomonas<-RelGenBU1[which(RelGenBU1$Genus == "Sphingomonas"),]; head(Sphingomonas$Genus)
Streptomyces<-RelGenBU1[which(RelGenBU1$Genus == "Streptomyces"),]; head(Streptomyces$Genus)


#Burned samples ......................................................................
Acidibacter$SampleType<-as.factor(Acidibacter$SampleType)
Bryobacter$SampleType<-as.factor(Bryobacter$SampleType)
Burkholderia$SampleType<-as.factor(Burkholderia$SampleType)
CanUndaeobacter$SampleType<-as.factor(CanUndaeobacter$SampleType)
Chitinophagaceae$SampleType<-as.factor(Chitinophagaceae$SampleType)
Chthoniobacter$SampleType<-as.factor(Chthoniobacter$SampleType)
Gemmatimonas$SampleType<-as.factor(Gemmatimonas$SampleType)
Massilia$SampleType<-as.factor(Massilia$SampleType)
Mucilaginibacter$SampleType<-as.factor(Mucilaginibacter$SampleType)
Mycobacterium$SampleType<-as.factor(Mycobacterium$SampleType)
Paenibacillus$SampleType<-as.factor(Paenibacillus$SampleType)
RB41$SampleType<-as.factor(RB41$SampleType)
Sphingomonas$SampleType<-as.factor(Sphingomonas$SampleType)
Streptomyces$SampleType<-as.factor(Streptomyces$SampleType)

Acidibacter$SampleType<-try(relevel(Acidibacter$SampleType,ref="Original"))
Bryobacter$SampleType<-try(relevel(Bryobacter$SampleType,ref="Original"))
Burkholderia$SampleType<-try(relevel(Burkholderia$SampleType,ref="Original"))
CanUndaeobacter$SampleType<-try(relevel(CanUndaeobacter$SampleType,ref="Original"))
Chitinophagaceae$SampleType<-try(relevel(Chitinophagaceae$SampleType,ref="Original"))
Chthoniobacter$SampleType<-try(relevel(Chthoniobacter$SampleType,ref="Original"))
Gemmatimonas$SampleType<-try(relevel(Gemmatimonas$SampleType,ref="Original"))
Massilia$SampleType<-try(relevel(Massilia$SampleType,ref="Original"))
Mycobacterium$SampleType<-try(relevel(Mycobacterium$SampleType,ref="Original"))
Mucilaginibacter$SampleType<-try(relevel(Mucilaginibacter$SampleType,ref="Original"))
Paenibacillus$SampleType<-try(relevel(Paenibacillus$SampleType,ref="Original"))
RB41$SampleType<-try(relevel(RB41$SampleType,ref="Original"))
Sphingomonas$SampleType<-try(relevel(Sphingomonas$SampleType,ref="Original"))
Streptomyces$SampleType<-try(relevel(Streptomyces$SampleType,ref="Original"))



#Explore Temporal correlations to ensure that we are capturing as much of the variance............... .................-----
#Unburned
ASVneg1 <- lmer(Abundance ~ SampleType + (1|Plot)+(1|TSFdays)+(1|TSFdays), data = Bryobacter)#model failed to converge
ASVneg2 <- lmer(Abundance ~ SampleType  +  (1|Plot)+(1|TSFdays), data = Bryobacter)
ASVneg3 <- lmer(Abundance ~ SampleType  + (1|Subplot)+(1|TSFdays), data = Bryobacter)
ASVneg4 <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Bryobacter)
ASVneg5 <- lmer(Abundance ~ SampleType  +  (1|Plot), data = Bryobacter)
ASVneg6 <- lmer(Abundance ~ SampleType  + (1|Subplot), data = Bryobacter)

#Check model to see which is better
AICc(ASVneg1,ASVneg2,ASVneg3,ASVneg4,ASVneg5,ASVneg6)#ASVneg4 better


#Begin statistical models.....................................................................................
AcidiLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Acidibacter)
Anova(AcidiLm, type=3)#0.3994 

BryoLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Bryobacter)
Anova(BryoLm, type=3)#0.008793 ** 
PHBryo<-summary(glht(BryoLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHBryo#T-O (0.022) 

BurkLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Burkholderia)
Anova(BurkLm, type=3)#0.310650 

CandiLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = CanUndaeobacter)
Anova(CandiLm, type=3)#0.697541 4 

ChitiLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Chitinophagaceae)
Anova(ChitiLm, type=3)#0.3994 
#PHChiti<-summary(glht(ChitiLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHChiti#T-O (0.000499) T-F(0.005632)

ChthoLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Chthoniobacter)
Anova(ChthoLm, type=3)#0.0002688 ***
PHChtho<-summary(glht(ChthoLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHChtho#T-O (0.000499) T-F(0.005632)

GemmaLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Gemmatimonas)
Anova(GemmaLm, type=3)#0.7545

MassiLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Massilia)
Anova(MassiLm, type=3)#0.3235

MuciLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Mucilaginibacter)
Anova(MuciLm, type=3)#0.0516381

MycoLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Mycobacterium)
Anova(MycoLm, type=3)#0.0129070 *
PHMyco<-summary(glht(MycoLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHMyco#T-O (0.0121)

PaeniLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Paenibacillus)
Anova(PaeniLm, type=3)#0.01503 *
PHPaeni<-summary(glht(PaeniLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHPaeni#T-O (0.0351);T-F (0.0351)

RB41Lm <- lmer(Abundance ~ SampleType  + (1|Plot), data = RB41)
Anova(RB41Lm, type=3)#0.2422956

SphinLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Sphingomonas)
Anova(SphinLm, type=3)#0.8643 

StrepLm <- lmer(Abundance ~ SampleType  + (1|Plot), data = Streptomyces)
Anova(StrepLm, type=3)#0.59083


