#Reset R's Brain
rm(list=ls())


#Set working directory................................................................................................
setwd("/Users/fabipc/Dropbox/2-UC-Riverside/1-Dissertation/1-Research/1-HolyFire/4-Freeze-thaw-Analysis/1-Lib7Thawed-NoCNF9WT2")


#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)
library(car)#anova for glmer

Metadata<-read_tsv("Fungi/Metadata/MetaRareFun-CNF9WT2-2.tsv")#Rare metadata
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


#--Subset data by treatment (burned vs unburned) ......................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) 

# * * Subset All to maintain only the burned samples of all data......................
physeqB<-subset_samples(physeq,Treatment=="Burned");physeqB#2607x27
physeqUn<-subset_samples(physeq,Treatment=="Unburned");physeqUn#2607x27



#.----
#********************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -----------------------------------------------------------------
#*********************************************************************************************-----

#* * Burned ........................................................................
GenFB<- tax_glom(physeqB, "Genus")
RelGenFB <- transform_sample_counts(GenFB, function(x) x / sum(x))
RelGenFB1<-psmelt(RelGenFB)

#* * Unburned ......................................................................
GenFU<- tax_glom(physeqUn, "Genus")
RelGenFU<- transform_sample_counts(GenFU, function(x) x / sum(x))
RelGenFU1<-psmelt(RelGenFU)




#.----
#.----
#**************************************************************************************************************************-----
# EXPORT DATA FILES --------------------------------------------------------------------------------------------------------------
#**************************************************************************************************************************-----
dir.create(file.path("Fungi/Analysis/RelativeAbundance/Tables/DeStats"), recursive=TRUE)


#Re-Import soil samples ......................................................................................................
write.csv(RelGenFB1,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/RelAbund-Burn-All.csv")
write.csv(RelGenFU1,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/RelAbund-Unburn-All.csv")


#Reimport soil samples
RelGenFB1<-read.csv("Fungi/Analysis/RelativeAbundance/Tables/DeStats/RelAbund-Burn-All.csv")
RelGenFU1<-read.csv("Fungi/Analysis/RelativeAbundance/Tables/DeStats/RelAbund-Unburn-All.csv")

#Quality control ..................................................................
RelGenFB1$SampleType<-as.factor(RelGenFB1$SampleType)
RelGenFU1$SampleType<-as.factor(RelGenFU1$SampleType)

RelGenFB1$SampleType<-try(relevel(RelGenFB1$SampleType,ref="Original"))
RelGenFU1$SampleType<-try(relevel(RelGenFU1$SampleType,ref="Original"))


#---
#-----
#********************************************************************************************************************----
#---DESCRIPTIVE STATISTICS -------------------------------------------------------------------------------------------------
#********************************************************************************************************************----
#Subset data only thawed and original 
FunBOT <- subset(RelGenFB1, SampleType == "Original" | SampleType == "Thawed");unique(FunBOT$SampleType)
FunBOF <- subset(RelGenFB1, SampleType == "Original" | SampleType == "Frozen");unique(FunBOF$SampleType)

FunUOT <- subset(RelGenFU1, SampleType == "Original" | SampleType == "Thawed");unique(FunUOT$SampleType)
FunUOF <- subset(RelGenFU1, SampleType == "Original" | SampleType == "Frozen");unique(FunUOF$SampleType)


#Relevel to make original the base level
FunBOT$SampleType<-try(relevel(FunBOT$SampleType,ref="Original"))
FunBOF$SampleType<-try(relevel(FunBOF$SampleType,ref="Original"))

FunUOT$SampleType<-try(relevel(FunUOT$SampleType,ref="Original"))
FunUOF$SampleType<-try(relevel(FunUOF$SampleType,ref="Original"))


#Burned orig-thaw
SampFBOT<-FunBOT %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, SampleType) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Mean = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));SampFBOT

#Orig-Frozen
SampFBOF<-FunBOF %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, SampleType) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Mean = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));SampFBOF



write.csv(SampFBOT,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-B-OrigThaw-DesStats-All.csv")
write.csv(SampFBOF,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-B-OrigFrozen-DesStats-All.csv")


#Same stats but for Unburned dataset
#Burned orig-thaw
SampFUOT<-FunUOT %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, SampleType) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Mean = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));SampFUOT

#Orig-Frozen
SampFUOF<-FunUOF %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, SampleType) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
              mutate(percent = round (percent, 2)));SampFUOF



write.csv(SampFUOT,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-U-OrigThaw-DesStats-All.csv")
write.csv(SampFUOF,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-U-OrigFrozen-DesStats-All.csv")


#----
#
#***********************************************************************************************************************----
#--- PERCENT CHANGE CALCULATIONS -------------------------------------------------------------------------------------------
#***********************************************************************************************************************----
library(dplyr)

#Burned orig-thaw
FBOT<-FunBOT %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();head(FBOT)

#Orig-Frozen
FBOF<-FunBOF %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();head(FBOF)



write.csv(FBOT,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-B-OrigThaw-PerChange-All.csv")
write.csv(FBOF,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-B-OrigFrozen-PerChange-All.csv")


#Same stats but for Unburned dataset
#Burned orig-thaw
FUOT<-FunUOT %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();FUOT

#Orig-Frozen
FUOF<-FunUOF %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,SampleType) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 3))%>%
  as.data.frame();FUOF



write.csv(FUOT,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-U-OrigThaw-PerChange-All.csv")
write.csv(FUOF,"Fungi/Analysis/RelativeAbundance/Tables/DeStats/Fun-U-OrigFrozen-PerChange-All.csv")





#----
#-----
#***********************************************************************************************************************-----
# CALCULATE SIGNIFICANCE----------------------------------------------------------------------------------------------------
#***********************************************************************************************************************-----
#Re-Import soil samples ......................................................................................................
RelGenFB1<-read.csv("Fungi/Analysis/RelativeAbundance/Tables/DeStats/RelAbund-Burn-All.csv")
RelGenFU1<-read.csv("Fungi/Analysis/RelativeAbundance/Tables/DeStats/RelAbund-Unburn-All.csv")



#SUBSET THE DATA FOR TREATMENT ........................................................
attach(RelGenFB1)

Coniochaeta<-RelGenFB1[which(RelGenFB1$Genus == "Coniochaeta"),]; head(Coniochaeta$Genus)
Penicillium<-RelGenFB1[which(RelGenFB1$Genus == "Penicillium"),]; head(Penicillium$Genus)
Pyronema<-RelGenFB1[which(RelGenFB1$Genus == "Pyronema"),]; head(Pyronema$Genus)
Fimetariella<-RelGenFB1[which(RelGenFB1$Genus == "Fimetariella"),]; head(Fimetariella$Genus)
Tomentella<-RelGenFB1[which(RelGenFB1$Genus == "Tomentella"),]; head(Tomentella$Genus)
Aspergillus<-RelGenFB1[which(RelGenFB1$Genus == "Aspergillus"),]; head(Aspergillus$Genus)
Geminibasidium<-RelGenFB1[which(RelGenFB1$Genus == "Geminibasidium"),]; head(Geminibasidium$Genus)
Coprinellus<-RelGenFB1[which(RelGenFB1$Genus == "Coprinellus"),]; head(Coprinellus$Genus)
Inocybe<-RelGenFB1[which(RelGenFB1$Genus == "Inocybe"),]; head(Inocybe$Genus)


#CHECK DATA FOR NORMALITY..........................................................................
"Burkholderia-Caballeronia-Paraburkholderia"

#Check histogram and run shapiro test ...........................................
hist(Coniochaeta$Abundance); shapiro.test(Coniochaeta$Abundance)
hist(Penicillium$Abundance); shapiro.test(Penicillium$Abundance)
hist(Pyronema$Abundance); shapiro.test(Pyronema$Abundance)
hist(Fimetariella$Abundance); shapiro.test(Fimetariella$Abundance)
hist(Tomentella$Abundance); shapiro.test(Tomentella$Abundance)

#All not normal...

#-----
#----
#*****************************************************************************************************************************************----
#SIGNIFICANCE PER TAXA PER TREATMENT----------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************----
library(MASS)
library(lmerTest)
library(lme4)
library(MuMIn)

#Quality control ..................................................................
Penicillium$SampleType<-as.factor(Penicillium$SampleType)
Penicillium$SampleType<-as.factor(Penicillium$SampleType)
Coniochaeta$SampleType<-as.factor(Coniochaeta$SampleType)
Fimetariella$SampleType<-as.factor(Fimetariella$SampleType)
Tomentella$SampleType<-as.factor(Tomentella$SampleType)
Aspergillus$SampleType<-as.factor(Aspergillus$SampleType)
Geminibasidium$SampleType<-as.factor(Geminibasidium$SampleType)
Coprinellus$SampleType<-as.factor(Coprinellus$SampleType)
Inocybe$SampleType<-as.factor(Inocybe$SampleType)


Aspergillus$SampleType<-try(relevel(Aspergillus$SampleType,ref="Original"))
Coniochaeta$SampleType<-try(relevel(Coniochaeta$SampleType,ref="Original"))
Pyronema$SampleType<-try(relevel(Pyronema$SampleType,ref="Original"))
Penicillium$SampleType<-try(relevel(Penicillium$SampleType,ref="Original"))
Fimetariella$SampleType<-try(relevel(Fimetariella$SampleType,ref="Original"))
Tomentella$SampleType<-try(relevel(Tomentella$SampleType,ref="Original"))
Geminibasidium$SampleType<-try(relevel(Geminibasidium$SampleType,ref="Original"))
Coprinellus$SampleType<-try(relevel(Coprinellus$SampleType,ref="Original"))
Inocybe$SampleType<-try(relevel(Inocybe$SampleType,ref="Original"))



#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
ASVneg1 <- lmer(Abundance ~ SampleType + (1|Plot)+(1|TSFdays)+(1|TSFdays), data = Pyronema)#model failed to converge
ASVneg2 <- lmer(Abundance ~ SampleType  +  (1|Plot)+(1|TSFdays), data = Pyronema)
ASVneg3 <- lmer(Abundance ~ SampleType  + (1|Subplot)+(1|TSFdays), data = Pyronema)
ASVneg4 <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Pyronema)
ASVneg5 <- lmer(Abundance ~ SampleType  +  (1|Plot), data = Pyronema)
ASVneg6 <- lmer(Abundance ~ SampleType  + (1|Subplot), data = Pyronema)

#Check model to see which is better
AICc(ASVneg1,ASVneg2,ASVneg3,ASVneg4,ASVneg5,ASVneg6)#ASVneg1 better
library(multcomp) #for post-hoc multilevel test

#Test for significance----------------------------------------------------------------------------------------------------
PyroLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Pyronema)
Anova(PyroLm, type=3)# sample type was signficant 0.002290 **
PHPyro<-summary(glht(PyroLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHPyro#T-O (0.00315) TF= 0.01517

PeniLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Penicillium)
Anova(PeniLm, type=3)# sample type was signficant 0.04982 *
PHPeni<-summary(glht(PeniLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHPeni#TF=0.0459 *

ConiLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Coniochaeta)
Anova(ConiLm, type=3)# sample type was signficant 0.30377 
PHConi<-summary(glht(ConiLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHConi#no signif

FimeLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Fimetariella)
Anova(FimeLm, type=3)# sample type was signficant 0.9366
PHFime<-summary(glht(FimeLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHFime#No signif

AsperLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Aspergillus)
Anova(AsperLm, type=3)#not significant 0.2140
PHAsper<-summary(glht(AsperLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHAsper#No signif

GeminLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Geminibasidium)
Anova(GeminLm, type=3)#not significant 0.4540
PHGemin<-summary(glht(GeminLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHGemin#No signif

CopriLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Coprinellus)
Anova(CopriLm, type=3)#not significant 0.3766
PHCopri<-summary(glht(CopriLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHCopri#No signif


TomenLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Tomentella)
Anova(TomenLm, type=3)#not significant 0.2638
PHTomen<-summary(glht(TomenLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHTomen#No signif

InocybeLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Inocybe)
Anova(InocybeLm, type=3)#not significant 0.1462
PHInocybe<-summary(glht(InocybeLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHInocybe#No signif




#-----
#----
#*****************************************************************************************************************************************----
#UNBURNED SAMPLES ----------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************----
attach(RelGenFU1)

Coniochaeta<-RelGenFU1[which(RelGenFU1$Genus == "Coniochaeta"),]; head(Coniochaeta$Genus)
Cortinarius<-RelGenFU1[which(RelGenFU1$Genus == "Cortinarius"),]; head(Cortinarius$Genus)
Cladophialophora<-RelGenFU1[which(RelGenFU1$Genus == "Cladophialophora"),]; head(Cladophialophora$Genus)
Cenococcum<-RelGenFU1[which(RelGenFU1$Genus == "Cenococcum"),]; head(Cenococcum$Genus)
Balsamia<-RelGenFU1[which(RelGenFU1$Genus == "Balsamia"),]; head(Balsamia$Genus)
Geopora<-RelGenFU1[which(RelGenFU1$Genus == "Geopora"),]; head(Geopora$Genus)
Hyaloscypha<-RelGenFU1[which(RelGenFU1$Genus == "Hyaloscypha"),]; head(Hyaloscypha$Genus)
Hyphodiscus<-RelGenFU1[which(RelGenFU1$Genus == "Hyphodiscus"),]; head(Hyphodiscus$Genus)
Inocybe<-RelGenFU1[which(RelGenFU1$Genus == "Inocybe"),]; head(Inocybe$Genus)
Penicillium<-RelGenFU1[which(RelGenFU1$Genus == "Penicillium"),]; head(Penicillium$Genus)
Penicillium<-RelGenFU1[which(RelGenFU1$Genus == "Penicillium"),]; head(Penicillium$Genus)
Ramaria<-RelGenFU1[which(RelGenFU1$Genus == "Ramaria"),]; head(Ramaria$Genus)
Thelephora<-RelGenFU1[which(RelGenFU1$Genus == "Thelephora"),]; head(Thelephora$Genus)
Tomentella<-RelGenFU1[which(RelGenFU1$Genus == "Tomentella"),]; head(Tomentella$Genus)
Trechispora<-RelGenFU1[which(RelGenFU1$Genus == "Trechispora"),]; head(Trechispora$Genus)
Venturia<-RelGenFU1[which(RelGenFU1$Genus == "Venturia"),]; head(Venturia$Genus)


#---
#---
#*****************************************************************************************************************************************----
#SIGNIFICANCE PER TAXA PER TREATMENT----------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************----


#Quality control ..................................................................
Cenococcum$SampleType<-as.factor(Cenococcum$SampleType)
Cortinarius$SampleType<-as.factor(Cortinarius$SampleType)
Cladophialophora$SampleType<-as.factor(Cladophialophora$SampleType)
Coniochaeta$SampleType<-as.factor(Coniochaeta$SampleType)
Balsamia$SampleType<-as.factor(Balsamia$SampleType)
Geopora$SampleType<-as.factor(Geopora$SampleType)
Hyaloscypha$SampleType<-as.factor(Hyaloscypha$SampleType)
Hyphodiscus$SampleType<-as.factor(Hyphodiscus$SampleType)
Inocybe$SampleType<-as.factor(Inocybe$SampleType)
Mycena$SampleType<-as.factor(Mycena$SampleType)
Penicillium$SampleType<-as.factor(Penicillium$SampleType)
Ramaria$SampleType<-as.factor(Ramaria$SampleType)
Thelephora$SampleType<-as.factor(Thelephora$SampleType)
Tomentella$SampleType<-as.factor(Tomentella$SampleType)
Trechispora$SampleType<-as.factor(Trechispora$SampleType)
Venturia$SampleType<-as.factor(Venturia$SampleType)


Balsamia$SampleType<-try(relevel(Balsamia$SampleType,ref="Original"))
Cenococcum$SampleType<-try(relevel(Cenococcum$SampleType,ref="Original"))
Cladophialophora$SampleType<-try(relevel(Cladophialophora$SampleType,ref="Original"))
Coniochaeta$SampleType<-try(relevel(Coniochaeta$SampleType,ref="Original"))
Cortinarius$SampleType<-try(relevel(Cortinarius$SampleType,ref="Original"))
Geopora$SampleType<-try(relevel(Geopora$SampleType,ref="Original"))
Hyaloscypha$SampleType<-try(relevel(Hyaloscypha$SampleType,ref="Original"))
Hyphodiscus$SampleType<-try(relevel(Hyphodiscus$SampleType,ref="Original"))
Inocybe$SampleType<-try(relevel(Inocybe$SampleType,ref="Original"))
Mycena$SampleType<-try(relevel(Mycena$SampleType,ref="Original"))
Penicillium$SampleType<-try(relevel(Penicillium$SampleType,ref="Original"))
Ramaria$SampleType<-try(relevel(Ramaria$SampleType,ref="Original"))
Thelephora$SampleType<-try(relevel(Thelephora$SampleType,ref="Original"))
Tomentella$SampleType<-try(relevel(Tomentella$SampleType,ref="Original"))
Trechispora$SampleType<-try(relevel(Trechispora$SampleType,ref="Original"))
Venturia$SampleType<-try(relevel(Venturia$SampleType,ref="Original"))





#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
ASVneg1 <- lmer(Abundance ~ SampleType + (1|Plot)+(1|TSFdays)+(1|TSFdays), data = Cenococcum)#model failed to converge
ASVneg2 <- lmer(Abundance ~ SampleType  +  (1|Plot)+(1|TSFdays), data = Cenococcum)
ASVneg3 <- lmer(Abundance ~ SampleType  + (1|Subplot)+(1|TSFdays), data = Cenococcum)
ASVneg4 <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Cenococcum)
ASVneg5 <- lmer(Abundance ~ SampleType  +  (1|Plot), data = Cenococcum)
ASVneg6 <- lmer(Abundance ~ SampleType  + (1|Subplot), data = Cenococcum)

#Check model to see which is better
AICc(ASVneg1,ASVneg2,ASVneg3,ASVneg4,ASVneg5,ASVneg6)#Mod 4 better
library(multcomp) #for post-hoc multilevel test

#Test for significance----------------------------------------------------------------------------------------------------
BalsamiaLm <- glmer(Abundance ~ SampleType  + (1|Subplot), data = Balsamia)
Anova(BalsamiaLm, type=3)# sample type was signficant 0.1077
PHBalsamia<-summary(glht(BalsamiaLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHBalsamia#T-O (0.00315) TF= 0.01517

CenococcumLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Cenococcum)
Anova(CenococcumLm, type=3)# sample type was signficant 5.093e-05
PHCenococcum<-summary(glht(CenococcumLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHCenococcum#TF=0.00274, T0 5.78e-05

CladoLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Cladophialophora)
Anova(CladoLm, type=3)# sample type was signficant 0.214917
PHClado<-summary(glht(CladoLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHClado#no signif

ConiLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Coniochaeta)
Anova(ConiLm, type=3)# sample type was signficant 0.33463 
PHConi<-summary(glht(ConiLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHConi#No signif

CortLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Cortinarius)
Anova(CortLm, type=3)#not significant  0.2997
PHCort<-summary(glht(CortLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHCort#No signif

GeoporaLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Geopora)
Anova(GeoporaLm, type=3)#not significant 0.6810
PHGeopora<-summary(glht(GeoporaLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHGeopora#No signif

HyaloLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Hyaloscypha)
Anova(HyaloLm, type=3)#not significant  0.09493 .
PHHyalo<-summary(glht(HyaloLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHHyalo#No signif

HyphoLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Hyphodiscus)
Anova(HyphoLm, type=3)#not significant 0.2991
PHHypho<-summary(glht(HyphoLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHHypho#No signif

InocybeLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Inocybe)
Anova(InocybeLm, type=3)#not significant 0.56022 
PHInocybe<-summary(glht(InocybeLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHInocybe#No signif

MyceLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Mycena)
Anova(MyceLm, type=3)#not significant 0.301310 
PHMyce<-summary(glht(MyceLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHMyce#No signif

PeniLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Penicillium)
Anova(PeniLm, type=3)#not significant  2.56e-08 ***
PHPeni<-summary(glht(PeniLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHPeni #TO = <0.00001; TF < 0.00001

RamaLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Ramaria)
Anova(RamaLm, type=3)#not significant 0.3253
PHRama<-summary(glht(RamaLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHRama #TO = <0.00001; TF < 0.00001

TheleLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Thelephora)
Anova(TheleLm, type=3)#not significant 0.1287
PHThele<-summary(glht(TheleLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHThele #TO = <0.00001; TF < 0.00001

TomenLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Tomentella)
Anova(TomenLm, type=3)#not significant 0.2067
PHTomen<-summary(glht(TomenLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHTomen #TO = <0.00001; TF < 0.00001

TrechiLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Trechispora)
Anova(TrechiLm, type=3)#not significant 0.2253
PHTrechi<-summary(glht(TrechiLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHTrechi #TO = <0.00001; TF < 0.00001

VenturiaLm <- lmer(Abundance ~ SampleType  + (1|Plot)+(1|TSFdays), data = Venturia)
Anova(VenturiaLm, type=3)#not significant 0.07996 .
PHVenturia<-summary(glht(VenturiaLm, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PHVenturia #TO = <0.00001; TF < 0.00001



# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt"), recursive=TRUE)

#Export ANOVA values ...............................................................................................................
capture.output(Anova(Bryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Bryobacterglmer.csv")
capture.output(Anova(Bryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Can.Udaeobacterglmer.csv")
capture.output(Anova(Udaeobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Sphingomonasglmer.csv")


---
#-----
#----------
#********************************************************************************************************************************----
# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/"), recursive = TRUE)

#Generalized Regression ....................................................................................................................
capture.output(Anova(AdhaeribacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Adhaeribacter-GlmerBurn.csv")
capture.output(Anova(BryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Bryobacter-GlmerBurn.csv")
capture.output(Anova(BlastococcusLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Blastococcus-GlmerBurn.csv")
capture.output(Anova(Can.UdaeobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Can.Udaeobacter-GlmerBurn.csv")
capture.output(Anova(RB41LM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/RB41-GlmerBurn.csv")
capture.output(Anova(ConexibacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Conexibacter-GlmerBurn.csv")
capture.output(Anova(DomiBryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/DomiBryobacter-GlmerBurn.csv")
capture.output(Anova(GemmatimonasLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Gemmatimonas-GlmerBurn.csv")
capture.output(Anova(MucilaginibacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Mucilaginibacter-GlmerBurn.csv")
capture.output(Anova(NoviherbaspirillumLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Noviherbaspirillum-GlmerBurn.csv")
capture.output(Anova(PaeniBryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/PaeniBryobacter-GlmerBurn.csv")
capture.output(Anova(PedobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Pedobacter-GlmerBurn.csv")
capture.output(Anova(RB41LM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/RB41-GlmerBurn.csv")
capture.output(Anova(SolirubrobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Solirubrobacter-GlmerBurn.csv")
capture.output(Anova(TumeBryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/TumeBryobacter-GlmerBurn.csv")
detach(RelGenSoilB1)


