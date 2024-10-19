#5-14-2022
 
rm(list=ls())#reset working directory
 
#Set working directory-----------------------------------------------------------------------------------
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/Lib7Thawed-NoCNF9WT2/Bacteria/")

#LOAD PACKAGES ------------------------------------------------------------------------------------
library("pairwiseAdonis") #test for significance, post-hoc 1-Analysis 
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)


# LOAD DATA (rarefied table and metadata)-------------------------------------------------------------------
MetaRareBac<- read.csv("Metadata/MetaRareSppBac-CNF9WT2.csv", na.strings = "NA", header = TRUE)
 
#.----
#.----
#**********************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#**********************************************************************************************************----
 
# * LOOK AT DATA STRUCTURE ..................................................
MetaRareBac$Plot<- as.factor(MetaRareBac$Plot)
MetaRareBac$Treatment<-as.factor(MetaRareBac$Treatment)
MetaRareBac$SampleType<-as.factor(MetaRareBac$SampleType)

#CHANGE BASE LEVEL ...........................................................  
MetaRareBac$Treatment <- try(relevel(MetaRareBac$Treatment , "Unburned"))
MetaRareBac$SampleType <- try(relevel(MetaRareBac$SampleType , "Original"))

#----
#----
#******************************************************************************************----
# CHECK DATA FOR NORMA----------------------------------------------------
#******************************************************************************************----
attach(MetaRareBac)
par(mfrow=c(1,1))
hist(S.obs)

#shapiro test 
shapiro.test(sqrt(MetaRareBac$S.obs))#0.344

#For statistical analysis--use the squareroot of the data
MetaRareBac$SqASV<-sqrt(MetaRareBac$S.obs)


#Determine percent change between treatment at different samples types..................
library(dplyr)
library(srvyr)#prevents error for dplyr

TrtSobs <-  MetaRareBac%>% 
  group_by(SampleType,Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtSobs

#-43.44289 (original), -37.29573(frozen), -62.94026 (thawed)
#Export ghaphs and data
write.csv(TrtSobs, "Analysis/Diversity/Tables/PercentChangeTrtSampleType-Bac.csv")




#----
#-----
#***********************************************************************************************************----
# ------------------ PLOTS SPECIES RICHNESS (ESV's)  ---------------------------------------
#***********************************************************************************************************----
#REATMENT -----------------------------------------------------------------
SppTrt<-ggplot(MetaRareBac, aes(x=Treatment, y=S.obs)) +
  geom_boxplot(aes(fill=Treatment), show.legend = F)+ 
  theme_bw() +  ylab("Species Richness") + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "bottom",
        text = element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16, angle = 90,  hjust=.5, colour = "black"),
        axis.text.x = element_text(size=16, color = c("#2f5e59","#664128")))+ ylim(0,1050)+
scale_fill_manual(values=c("#45877f","#a2673f"));SppTrt

SppTrt2<-SppTrt + facet_wrap(~SampleType)+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=14));SppTrt2

SppTrt3<-SppTrt2 + stat_compare_means(method = "anova", size=5,label.x=1.45,label.y = 1035)+
  annotate("text", x=2.3,  y=1000,label = "-43%", size=5); SppTrt3


#Export the graphs................................................................................
dir.create(file.path("Analysis/Diversity/Graphs"), recursive = TRUE)

pdf("Analysis/Diversity/Graphs/Alpha-Bac-Trt-SampleType-FOT-CNF9WT2.pdf",height=5, width=9,onefile=FALSE)
SppTrt3
dev.off()





#---------------------------------------------------------------------------------------------------------
#Look at species richness only comparing the bunred and unburned samples independently
SppIndBac<-ggplot(MetaRareBac, aes(x=SampleType, y=S.obs)) +
  geom_boxplot(aes(fill=SampleType),show.legend = F)+ 
  theme_bw() +  ylab("Species Richness") + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "bottom",
        text = element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=16, angle = 90, hjust=0.5,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        legend.text = element_text(size=18))+ ylim(0,1350)+
  scale_fill_manual(values=c("#d8dbed","#8c91a8","#4c5478"));SppIndBac

SppIndBac2<-SppIndBac + facet_wrap(~Treatment)+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=14));SppIndBac2

my_comparisons <- list(c("Original", "Thawed"),c("Original", "Frozen"),c("Thawed", "Frozen"))
SppIndBac3<-SppIndBac2 + stat_compare_means(comparisons = my_comparisons);SppIndBac3

-
pdf("Analysis/Diversity/Graphs/Alpha-Bac-Total-Burned-Unburned-FOT-CNF9WT2.pdf",height=5, width=8,onefile=FALSE)
SppIndBac3
dev.off()


#Set original as the base level
#Subset data and test for significance-----------------------------------------
attach(MetaRareBac)
BurBac<-MetaRareBac[which(Treatment=="Burned"),];head(BurBac[10:16])
UnBac<-MetaRareBac[which(Treatment=="Unburned"),];head(UnBac[10:16])


STSobsB <- BurBac%>% 
  group_by(SampleType) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();STSobsB

STSobsUn <- UnBac%>% 
  group_by(SampleType) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();STSobsUn


#-43.44289 (original), -37.29573(frozen), -62.94026 (thawed)
#Export ghaphs and data
write.csv(STSobsB, "Analysis/Diversity/Tables/PercentChangeSampleType-Bac-Burned.csv")
write.csv(STSobsUn, "Analysis/Diversity/Tables/PercentChangeSampleType-Bac-Unburned.csv")





#----
#----
#*************************************************************************************----
#Test significance
#*************************************************************************************----
library(multcomp) #for post-hoc multilevel test

#Test nested level
attach(BurBac)
Mod1<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot)+(1|TSFdays), data=BurBac);summary(Mod1)
Mod2<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot), data=BurBac);summary(Mod2)
Mod3<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=BurBac);summary(Mod3)
Mod4<-lmer(S.obs~SampleType+(1|Plot), data=BurBac);summary(Mod4)
Mod5<-lmer(S.obs~SampleType+(1|Subplot), data=BurBac);summary(Mod5)

#Check at models to see which one is better
anova(Mod1,Mod2,Mod3,Mod4,Mod5) #Mod 3, lower AIC

#Results
AovBBac<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=BurBac);summary(AovBBac)#no signif
PostHocB<-summary(glht(AovBBac, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PostHocB
detach(BurBac)

#Test nested level for the unburned samples.................................................
attach(UnBac)
Mod1U<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot)+(1|TSFdays), data=UnBac);summary(Mod1U)
Mod2U<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot), data=UnBac);summary(Mod2U)
Mod3U<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=UnBac);summary(Mod3U)
Mod4U<-lmer(S.obs~SampleType+(1|Plot), data=UnBac);summary(Mod4U)
Mod5U<-lmer(S.obs~SampleType+(1|Subplot), data=UnBac);summary(Mod5U)

#Check at models to see which one is better
anova(Mod1U,Mod2U,Mod3U,Mod4U,Mod5U) #Mod4 slightly better, but results sames as mod3 so for 
#consistency we are stickig to Mod3

#Results-------------------------------------------------------------------------------------------------------------------------------
AovUnBac<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=UnBac);summary(AovUnBac)#NoSignif
PostHocUn<-summary(glht(AovUnBac, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PostHocUn#thawed-Frozen


#Export avo restults
capture.output(summary(AovBBac),file="Analysis/Diversity/Tables/Alpha-Burned-Aov-Bac.csv")
capture.output(summary(AovUnBac),file="Analysis/Diversity/Tables/Alpha-Unburned-Aov-Bac.csv")
capture.output(PostHocB,file="Analysis/Diversity/Tables/Alpha-PostHoc-Burned-Bac.csv")
capture.output(PostHocUn,file="Analysis/Diversity/Tables/Alpha-PostHoc-Unburned-Bac.csv")



#----
#----
#*************************************************************************************----
#Treatment effect
#*************************************************************************************----
attach(MetaRareBac)
Orig<-MetaRareBac[which(SampleType=="Original"),];head(Orig[10:16])
Thaw<-MetaRareBac[which(SampleType=="Thawed"),];head(Orig[10:16])
Frozen<-MetaRareBac[which(SampleType=="Frozen"),];head(Orig[10:16])




#Test nested level
attach(Orig)
Mod1<-lmer(S.obs~Treatment+(1|Plot)+(1|Subplot)+(1|TSFdays), data=Orig);summary(Mod1)
Mod2<-lmer(S.obs~Treatment+(1|Plot)+(1|Subplot), data=Orig);summary(Mod2)
Mod3<-lmer(S.obs~Treatment+(1|Plot)+(1|TSFdays), data=Orig);summary(Mod3)
Mod4<-lmer(S.obs~Treatment+(1|Plot), data=Orig);summary(Mod4)
Mod5<-lmer(S.obs~Treatment+(1|Subplot), data=Orig);summary(Mod5)

#Thawed
Mod1T<-lmer(S.obs~Treatment+(1|Plot)+(1|Subplot)+(1|TSFdays), data=Thaw);summary(Mod1T)
Mod2T<-lmer(S.obs~Treatment+(1|Plot)+(1|Subplot), data=Thaw);summary(Mod2T)
Mod3T<-lmer(S.obs~Treatment+(1|Plot)+(1|TSFdays), data=Thaw);summary(Mod3T)
Mod4T<-lmer(S.obs~Treatment+(1|Plot), data=Thaw);summary(Mod4T)
Mod5T<-lmer(S.obs~Treatment+(1|Subplot), data=Thaw);summary(Mod5T)

#Frozen
Mod1F<-lmer(S.obs~Treatment+(1|Plot)+(1|Subplot)+(1|TSFdays), data=Frozen);summary(Mod1F)
Mod2F<-lmer(S.obs~Treatment+(1|Plot)+(1|Subplot), data=Frozen);summary(Mod2F)
Mod3F<-lmer(S.obs~Treatment+(1|Plot)+(1|TSFdays), data=Frozen);summary(Mod3F)
Mod4F<-lmer(S.obs~Treatment+(1|Plot), data=Frozen);summary(Mod4F)
Mod5F<-lmer(S.obs~Treatment+(1|Subplot), data=Frozen);summary(Mod5F)



#Check at models to see which one is better
anova(Mod1,Mod2,Mod3,Mod4,Mod5) #Mod 4, lower AIC
anova(Mod1T,Mod2T,Mod3T,Mod4T,Mod5T) #Mod 4, lower AIC
anova(Mod1F,Mod2F,Mod3F,Mod4F,Mod5F) #Mod 4, lower AIC


#Results Original
AovBacOrig<-lmer(S.obs~Treatment+(1|Plot), data=Orig);summary(AovBacOrig)#0.0951
PostHocOrig<-summary(glht(AovBacOrig, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"));PostHocOrig

#Results Thawed
AovBacThaw<-lmer(S.obs~Treatment+(1|Plot), data=Thaw);summary(AovBacThaw)#0.000829 ***
PostHocThaw<-summary(glht(AovBacThaw, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"));PostHocThaw

#Results Frozen
AovBacFrozen<-lmer(S.obs~Treatment+(1|Plot), data=Frozen);summary(AovBacFrozen)#0.0593 .
PostHocFrozen<-summary(glht(AovBacFrozen, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"));PostHocFrozen


#Export avo restults
capture.output(summary(AovBacOrig),file="Analysis/Diversity/Tables/Alpha-Trt-Aov-Orig-Bac.csv")
capture.output(summary(AovBacThaw),file="Analysis/Diversity/Tables/Alpha-Trt-Aov-Thaw-Bac.csv")
capture.output(summary(AovBacFrozen),file="Analysis/Diversity/Tables/Alpha-Trt-Aov-Frozen-Bac.csv")










