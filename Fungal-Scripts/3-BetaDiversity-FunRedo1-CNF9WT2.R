
rm(list=ls())#reset working directory
setwd("~/Dropbox/2-UC-Riverside/1-Dissertation/1-Research/1-HolyFire/4-Freeze-thaw-Analysis/1-Lib7Thawed-NoCNF9WT2/Fungi/")


#LOAD IBRARY ............................................................................................
library(ggplot2)
library(vegan)
library(pairwiseAdonis)# post-hoc analysis (or contrast) for adonis


#Load data..............................................................................................
MetaRareFun<- read.csv("Metadata/MetaRareFun-CNF9WT2.csv", na.strings = "N/A", row.names = 1,header = TRUE)
OtuTrans2<-read.csv("ASVtables/OtuTrans2-CNF9WT2.csv",  row.names = 1, check.names = FALSE)
OtuRareNZ<-read.csv("ASVtables/OtuRare3NZ-CNF9WT2.csv", row.names = 1, check.names = FALSE)


#----
#----
#*********************************************************************************************************----
#.........................QUALITY CONTROL.....................................................................
#*********************************************************************************************************----
all(rownames(MetaRareFun)== rownames(OtuTrans2))#makes sure names match

#DATA STRUCTURE .....................................
attach(MetaRareFun)
MetaRareFun$Plot <- as.factor(MetaRareFun$Plot)
MetaRareFun$SampleName <- as.factor(MetaRareFun$SampleName)
detach(MetaRareFun)


#Verify that names match
rownames(OtuRareNZ)[1]==rownames(MetaRareFun)[1]
dim(OtuRareNZ);dim(MetaRareFun)



#.----
#.----
#****************************************************************************************************----
# ------------------------BETA DIVERSITY (NMDS)----------------------------------------------------------
#****************************************************************************************************----

#All distance ...................................................................................
Dist<-avgdist(OtuTrans2, 5786, iterations=750, meanBac=median,transf= sqrt, dmethod="bray")
NMDS<-metaMDS(Dist, autotransform = FALSE, engine = "monoMDS", k=3, 
      weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDS)#Stress=0.0723742 




# * * Goodnes of fit-test .................................................
par(mfrow=c(1,1))
gof1<-goodness(NMDS);gof1

#Plot data to see quickly where points land................................
plot(NMDS, display="sites",cex=1.3)
stressplot(NMDS, p.col= "blue", l.col="red", lwd=2)#0.994, 0.968



#----
#----
#*******************************************************************************************************----
# .................NMDS PLOTS ............................................................................
#*******************************************************************************************************----

# * Check that names match in NMDS and metadata

#Extract NMDS scores (x and y coordinates)......................................
DataScores <- as.data.frame(scores(NMDS))


# * * EXPORT SCORES FOR PLOTTING DATA ..........................................
dir.create(file.path("Analysis/Diversity/Tables"), recursive = TRUE)
write.csv(DataScores, "Analysis/Diversity/Tables/DataScoresAll-CNF9WT2.csv")


# * * * To run graphs only without running distance again ......................
#Reimport datascore files for "all" and "soil only" ............................
DataScores<-read.csv("Analysis/Diversity/Tables/DataScoresAll-CNF9WT2.csv")

# * * * * Add columns to data frame to plot ....................................
DataScores$Treatment<- MetaRareFun$Treatment
DataScores$SampleType<- MetaRareFun$SampleType
DataScores$Plot<- MetaRareFun$Plot
DataScores$SampleName<- MetaRareFun$SampleName
DataScores$Timepoint<- MetaRareFun$Timepoint
DataScores$SampleID<- MetaRareFun$SampleID


#----
#----
#****************************************************************************************************----
#-------------------------------- CREATE GRAPHS ---------------------------------------------------------
#****************************************************************************************************----        
attach(MetaRareFun)
#col<-c("#d3e7f2","#80c5e8","#2a5bbd","#c4d4b8","#678764","#084d08",
#       "#f0d5ec","#cc93c3","#ab2494","#fff6c4","#f7e997","#ffd500",
#       "#a17777","#ba1e1e","#610505","#e3dede","#b8b8b8","#e8d7c8",
#       "#dea36f","#1a1918")


col<-c("#ede85f","#d3e6f0","#80c5e8","#32c7cf","#113e99",
       "#cde8ca","#8ba688","#308247","#ebe2be",
       "#ffae00","#856c3c","#850001","#d41111",
       "#c97171","#e3c5c5","#c2bcbc","#857f7f", 
       "#d9c1e6","#7309b5","#c97d3a")

#Treatment---thawed................................................................
attach(DataScores)
NMDSTrtSite<-ggplot(DataScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(shape = as.factor(SampleType), colour = SampleID), show.legend = T)+ 
  scale_shape_manual(values=c(19,17,15))+
  #scale_colour_manual(values = col)+
  theme(panel.background = element_blank(),
        legend.key.size = unit(5, "mm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size =12,colour = "black", angle=90,hjust=.5), 
        axis.text.x = element_text(size =12,colour = "black"), 
        axis.title = element_text(size = 18, colour = "black", face="bold"),
        legend.key=element_blank(),legend.position = "right",
        legend.title = element_text(size = 12, colour = "black",face = "bold"), 
        legend.text = element_text(size = 10, colour ="black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5))+
  guides(shape = guide_legend(override.aes = list(size = 4)))+ylim(-0.5,0.5)+
  labs(x = "NMDS1", colour = "Sample Name", y = "NMDS2", shape = "Sample Type");NMDSTrtSite

NMDSTrtSite2<-NMDSTrtSite + facet_wrap(~SampleID, )+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=13, face = "bold"));NMDSTrtSite2



#Treatment burned vs unburned.........................................................
NMDSTrt<-ggplot(DataScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(shape = as.factor(SampleType), colour = Treatment))+ 
  scale_shape_manual(values=c(19,17,15))+
  scale_colour_manual(values = c("#664027","#45877f"))+
  theme(panel.background = element_blank(),
        legend.key.size = unit(5, "mm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size =14,colour = "black", angle=90,hjust=.5), 
        axis.text.x = element_text(size =14,colour = "black"), 
        axis.title = element_text(size = 18, colour = "black", face="bold"),
        legend.key=element_blank(),legend.position = "right",
        legend.title = element_text(size = 14, colour = "black",face = "bold"), 
        legend.text = element_text(size = 12, colour ="black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5))+
  guides(shape = guide_legend(override.aes = list(size = 6)))+ ylim(-0.5,0.5)+
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2", shape = "Sample Type");NMDSTrt

NMDSTrt2<-NMDSTrt + facet_wrap(~Timepoint, scales = "free_y")+
  theme(strip.background =element_blank(),
          strip.text = element_text(size=14, face = "bold"));NMDSTrt2

#Just curious which points are mixing in T2
names<-plot(NMDS, type = "t")




#EXPORT GRAPHS ..........................................................................
pdf("Analysis/Diversity/Graphs/NMDS-SampleName-OFT-Fun-CNF9WT2.pdf", width=8, height=6)
NMDSTrtSite
dev.off()

pdf("Analysis/Diversity/Graphs/NMDS-SampleName-OFT-FacetFun-CNF9WT2-legend.pdf", width=14, height=10)
NMDSTrtSite2
dev.off()


#Treatment effects & facewrap info...........................................................
pdf("Analysis/Diversity/Graphs/NMDS-Sites-TrT-OFT-Fun-CNF9WT2.pdf", width=8, height=6)
NMDSTrt
dev.off()

pdf("Analysis/Diversity/Graphs/NMDS-TrT-facetTP-OFT-Fun-CNF9WT2.pdf", width=7.5, height=5)
NMDSTrt2
dev.off()




#----
#----
#**********************************************************************************----
#- .....Calculate significance ........................................................
##*********************************************************************************----
attach(MetaRareFun)

Adon<-adonis2(OtuRareNZ~ SampleType, method = "bray",permutations = 9999,
              data=MetaRareFun);Adon # R2= 0.01316, P= 0.36 

#Do contrast only is adonis is significant..............................................
#AdonPW<-pairwise.adonis2(OtuRareNZ ~Treatment2, data = MetaRareFun);AdonPW
capture.output(Adon,file="Analysis/Diversity/Tables/Adonis-Thawed-Original-Fun-CNF9WT2.csv")

#----
#----
#****************************************************************************************************----
#........Calculate significance per sample type ........................................................
#***************************************************************************************************----
#Import files 
OtuRareNZO<-read.csv("ASVtables/OtuRare3NZ-Orig.csv", row.names = 1, check.names = FALSE)
OtuRareNZT<-read.csv("ASVtables/OtuRare3NZ-Thaw.csv", row.names = 1, check.names = FALSE)
OtuRareNZF<-read.csv("ASVtables/OtuRare3NZ-Frozen.csv", row.names = 1, check.names = FALSE)


#Subset Metadata by samples type..................................................................
MetaO<-subset(MetaRareFun,SampleType=="Original");head(MetaO[,1:5])
MetaT<-subset(MetaRareFun,SampleType=="Thawed");head(MetaT[,1:5])
MetaF<-subset(MetaRareFun,SampleType=="Frozen");head(MetaF[,1:5])

#Statistical significance.........................................................................
AdonO<-adonis2(OtuRareNZO~ Treatment, method = "bray",permutations = 9999,
               data=MetaO);AdonO #R2=0.16329, P= 1e-04

AdonT<-adonis2(OtuRareNZT~ Treatment, method = "bray",permutations = 9999,
               data=MetaT);AdonT #R2=0.18476, P=2e-04 

AdonF<-adonis2(OtuRareNZF~ Treatment, method = "bray",permutations = 9999,
               data=MetaF);AdonF#R2=0.13373, P=3e-04  



