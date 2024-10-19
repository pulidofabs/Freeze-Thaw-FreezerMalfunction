#5-12-2022
 
rm(list=ls())#reset working directory
 
#Set working directory-----------------------------------------------------------------------------------
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/Lib7Thawed-NoCNF9WT2/Fungi")

#LOAD PACKAGES ------------------------------------------------------------------------------------
library("pairwiseAdonis") #test for significance, post-hoc 1-Analysis 
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)


# LOAD DATA (rarefied table and metadata)-------------------------------------------------------------------
MetaRareFun<- read.csv("Metadata/MetaRareSppFun-CNF9WT2.csv", na.strings = "NA", header = TRUE)
 
#.----
#.----
#**********************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#**********************************************************************************************************----
 
# * LOOK AT DATA STRUCTURE ..................................................
MetaRareFun$Plot<- as.factor(MetaRareFun$Plot)
MetaRareFun$Treatment<-as.factor(MetaRareFun$Treatment)
MetaRareFun$SampleType<-as.factor(MetaRareFun$SampleType)

#CHANGE BASE LEVEL ...........................................................  
MetaRareFun$Treatment <- try(relevel(MetaRareFun$Treatment , "Unburned"))
MetaRareFun$SampleType <- try(relevel(MetaRareFun$SampleType , "Original"))

#----
#----
#******************************************************************************************----
# CHECK DATA FOR NORMA----------------------------------------------------
#******************************************************************************************----
attach(MetaRareFun)
par(mfrow=c(1,1))
hist(S.obs)

#shapiro test 
shapiro.test(sqrt(MetaRareFun$S.obs))#0.3418 normal will do anove


#Determine percent change between treatment at different samples types..................
library(dplyr)
library(srvyr)#prevents error for dplyr

TrtSobs <-  MetaRareFun%>% 
  group_by(SampleType,Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtSobs

#-72.25111 (original), -66.46948 (frozen), -71.09609 (thawed)
#Export ghaphs and data
#write.csv(TrtSobs, "Analysis/Diversity/Tables/PercentChangeTrtSampleType-Fun.csv")


#----
#-----
#***********************************************************************************************************----
# ------------------ PLOTS SPECIES RICHNESS (ESV's)  ---------------------------------------
#***********************************************************************************************************----
#REATMENT -----------------------------------------------------------------
SppTrt<-ggplot(MetaRareFun, aes(x=Treatment, y=S.obs)) +
  geom_boxplot(aes(fill=Treatment), show.legend = F)+ 
  theme_bw() +  ylab("Species Richness") + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "bottom",
        text = element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16, angle = 90,  hjust=.9, colour = "black"),
        axis.text.x = element_text(size=16, color = c("#2f5e59","#664128")))+ 
  ylim(0,600)+scale_fill_manual(values=c("#45877f","#a2673f"));SppTrt

SppTrt2<-SppTrt + facet_wrap(~SampleType)+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=14));SppTrt2

SppTrt3<-SppTrt2 + stat_compare_means(method = "anova", size=5, 
              label.x=1.45,label.y = 590)+
  annotate("text", x=2.4,  y=569,label = "-72%", size=5); SppTrt3


#Export the graphs................................................................................
#dir.create(file.path("Analysis/Diversity/Graphs"), recursive = TRUE)

#pdf("Analysis/Diversity/Graphs/Alpha-Fun-Trt-SampleType-FOT-CNF9WT2.pdf",height=5, width=9,onefile=FALSE)
#SppTrt3
#dev.off()





#---------------------------------------------------------------------------------------------------------
#Look at species richness only comparing the bunred and unburned samples independently
SppIndFun<-ggplot(MetaRareFun, aes(x=SampleType, y=S.obs)) +
  geom_boxplot(aes(fill=SampleType))+ 
  theme_bw() +  ylab("Species Richness") + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "bottom",
        text = element_text(size=16),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=16, angle = 90, hjust=0.5,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"),
        legend.text = element_text(size=18))+ ylim(0,700)+
  scale_fill_manual(values=c("#d8dbed","#8c91a8","#4c5478"));SppIndFun

SppIndFun2<-SppIndFun + facet_wrap(~Treatment)+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=14));SppIndFun2

my_comparisons <- list(c("Original", "Thawed"),c("Original", "Frozen"),c("Thawed", "Frozen"))
SppIndFun3<-SppIndFun2 + stat_compare_means(comparisons = my_comparisons);SppIndFun3


#pdf("Analysis/Diversity/Graphs/Alpha-Fun-Total-Burned-Unburned-FOT-CNF9WT2.pdf",height=5, width=8,onefile=FALSE)
#SppIndFun3
#dev.off()


#Set original as the base level
#Subset data and test for significance-----------------------------------------
attach(MetaRareFun)
BurFun<-MetaRareFun[which(Treatment=="Burned"),];head(BurFun[10:16])
UnFun<-MetaRareFun[which(Treatment=="Unburned"),];head(UnFun[10:16])



SampObsB<-  BurFun%>% 
  group_by(SampleType) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();SampObsB



SampObsUn<-  UnFun%>% 
  group_by(SampleType) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();SampObsUn
#----
#----
#*************************************************************************************----
#Test significance
#*************************************************************************************----
library(multcomp) #for post-hoc multilevel test

#Test nested level
attach(BurFun)
Mod1<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot)+(1|TSFdays), data=BurFun);summary(Mod1)
Mod2<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot), data=BurFun);summary(Mod2)
Mod3<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=BurFun);summary(Mod3)
Mod4<-lmer(S.obs~SampleType+(1|Plot), data=BurFun);summary(Mod4)
Mod5<-lmer(S.obs~SampleType+(1|Subplot), data=BurFun);summary(Mod5)

#Check at models to see which one is better
anova(Mod1,Mod2,Mod3,Mod4,Mod5) #Mod 3, lower AIC

#Results
AovBFun<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=BurFun);summary(AovBFun)
       #Frozen different than original (P=0.00170)
PostHocB<-summary(glht(AovBFun, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PostHocB


#Test nested level for the unburned samples.................................................
attach(BurFun)
Mod1U<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot)+(1|TSFdays), data=UnFun);summary(Mod1U)
Mod2U<-lmer(S.obs~SampleType+(1|Plot)+(1|Subplot), data=UnFun);summary(Mod2U)
Mod3U<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=UnFun);summary(Mod3U)
Mod4U<-lmer(S.obs~SampleType+(1|Plot), data=UnFun);summary(Mod4U)
Mod5U<-lmer(S.obs~SampleType+(1|Subplot), data=UnFun);summary(Mod5U)

#Check at models to see which one is better
anova(Mod1U,Mod2U,Mod3U,Mod4U,Mod5U) #Mod4 better, but results sames as mod3 so for 
#consistency we are stickig to Mod3

#Results
AovUnFun<-lmer(S.obs~SampleType+(1|Plot)+(1|TSFdays), data=UnFun);summary(AovUnFun)#Frozen different than original (P=0.00599)
PostHocUn<-summary(glht(AovUnFun, linfct = mcp(SampleType = "Tukey")), test = adjusted("holm"));PostHocUn


#Export avo restults
capture.output(summary(AovBFun),file="Analysis/Diversity/Tables/Alpha-Burned-Aov-Fun.csv")
capture.output(summary(AovUnFun),file="Analysis/Diversity/Tables/Alpha-Unburned-Aov-Fun.csv")
capture.output(PostHocB,file="Analysis/Diversity/Tables/Alpha-PostHoc-Burned-Fun.csv")
capture.output(PostHocUn,file="Analysis/Diversity/Tables/Alpha-PostHoc-Unburned-Fun.csv")



#----
#----
#*************************************************************************************----
#Treatment effect
#*************************************************************************************----
attach(MetaRareFun)
Orig<-MetaRareFun[which(SampleType=="Original"),];head(Orig[10:16])
Thaw<-MetaRareFun[which(SampleType=="Thawed"),];head(Orig[10:16])
Frozen<-MetaRareFun[which(SampleType=="Frozen"),];head(Orig[10:16])




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
anova(Mod1,Mod2,Mod3,Mod4,Mod5) #Mod 3, lower AIC
anova(Mod1T,Mod2T,Mod3T,Mod4T,Mod5T) #Mod 3, lower AIC
anova(Mod1F,Mod2F,Mod3F,Mod4F,Mod5F) #Mod 3, lower AIC





#Results Original
AovFunOrig<-lmer(S.obs~Treatment+(1|Plot)+(1|TSFdays), data=Orig);summary(AovFunOrig)#(P=9.00e-05)
PostHocOrig<-summary(glht(AovFunOrig, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"));PostHocOrig

#Results Thawed
AovFunThaw<-lmer(S.obs~Treatment+(1|Plot)+(1|TSFdays), data=Thaw);summary(AovFunThaw)#(P=9.00e-05)
PostHocThaw<-summary(glht(AovFunThaw, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"));PostHocThaw

#Results Frozen
AovFunFrozen<-lmer(S.obs~Treatment+(1|Plot)+(1|TSFdays), data=Frozen);summary(AovFunFrozen)#(P=9.00e-05)
PostHocFrozen<-summary(glht(AovFunFrozen, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"));PostHocFrozen


#Export avo restults
capture.output(summary(AovFunOrig),file="Analysis/Diversity/Tables/Alpha-Trt-Aov-Orig-Fun.csv")
capture.output(summary(AovFunThaw),file="Analysis/Diversity/Tables/Alpha-Trt-Aov-Thaw-Fun.csv")
capture.output(summary(AovFunFrozen),file="Analysis/Diversity/Tables/Alpha-Trt-Aov-Frozen-Fun.csv")










