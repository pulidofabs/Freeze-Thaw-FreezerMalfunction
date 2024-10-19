rm(list=ls())#reset working directory

#Set working directory-----------------------------------------------------------------------------------
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/2-Freeze-thaw-Analysis/Lib7Thawed-NoCNF9WT2/Bacteria/")

#LOAD IBRARY ............................................................................................
library(ggplot2)
library(vegan)
library(pairwiseAdonis)# post-hoc analysis (or contrast) for adonis


#Load data..............................................................................................
MetaRareBac<- read.csv("Metadata/MetaRareBac-CNF9WT2.csv", na.strings = "N/A", row.names = 1,header = TRUE)
OtuTrans2<-read.csv("ASVtables/OtuTrans2-CNF9WT2.csv",  row.names = 1, check.names = FALSE)
OtuRareNZ<-read.csv("ASVtables/OtuRare3NZ-CNF9WT2.csv", row.names = 1, check.names = FALSE)


#----
#----
#*********************************************************************************************************----
#.........................QUALITY CONTROL.....................................................................
#*********************************************************************************************************----
all(rownames(MetaRareBac)==rownames(OtuTrans2))#makes sure names match

#DATA STRUCTURE .....................................
attach(MetaRareBac)
MetaRareBac$Plot <- as.factor(MetaRareBac$Plot)
MetaRareBac$SampleName <- as.factor(MetaRareBac$SampleName)
detach(MetaRareBac)


#Verify that names match
rownames(OtuRareNZ)[1]==rownames(MetaRareBac)[1]
dim(OtuRareNZ);dim(MetaRareBac)



#.----
#.----
#****************************************************************************************************----
# ------------------------BETA DIVERSITY (NMDS)----------------------------------------------------------
#****************************************************************************************************----

#All distance ...................................................................................
Dist<-avgdist(OtuTrans2, 6493, iterations=750, meanBac=median,transf= sqrt, dmethod="bray")
NMDS<-metaMDS(Dist, autotransform = FALSE, engine = "monoMDS", k=3, 
      weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDS)#Stress=0.07424754 




# * * Goodnes of fit-test .................................................
par(mfrow=c(1,1))
gof1<-goodness(NMDS);gof1

#Plot data to see quickly where points land................................
plot(NMDS, display="sites",cex=1.3)
stressplot(NMDS, p.col= "blue", l.col="red", lwd=2)#0.994, 0.976



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
DataScores$Treatment<- MetaRareBac$Treatment
DataScores$SampleType<- MetaRareBac$SampleType
DataScores$Plot<- MetaRareBac$Plot
DataScores$SampleName<- MetaRareBac$SampleName
DataScores$Timepoint<- MetaRareBac$Timepoint
DataScores$SampleID<- MetaRareBac$SampleID


#----
#----
#****************************************************************************************************----
#-------------------------------- CREATE GRAPHS ---------------------------------------------------------
#****************************************************************************************************----        
attach(MetaRareBac)
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
NMDSTrtSite<-ggplot(DataScores, aes(x = NMDS1, y = NMDS2, colour = SampleID)) + 
  geom_point(size = 4, aes(shape = as.factor(SampleType)), show.legend = T)+ 
  scale_shape_manual(values=c(19,17,15))+
  #scale_colour_manual(values = col)+
  theme(panel.background = element_blank(),
        legend.key.size = unit(5, "mm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size =11,colour = "black", angle=90,hjust=.5), 
        axis.text.x = element_text(size =11,colour = "black"), 
        axis.title = element_text(size = 16, colour = "black", face="bold"),
        legend.key=element_blank(),legend.position = "right",
        legend.title = element_text(size = 12, colour = "black",face = "bold"), 
        legend.text = element_text(size = 10, colour ="black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5))+
  guides(shape = guide_legend(override.aes = list(size = 4)))+ylim(-0.5,0.5)+
  labs(x = "NMDS1", y = "NMDS2", shape = "Sample Type");NMDSTrtSite

NMDSTrtSite2<-NMDSTrtSite + facet_wrap(~SampleID)+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=12, face = "bold"));NMDSTrtSite2

#FACET-WRAP BY SAMPLE TYPE
NMDSTrtSite<-ggplot(DataScores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = Treatment), show.legend = F)+ 
  scale_colour_manual(values = c("#664027","#45877f"))+
   scale_shape_manual(values=c(19,17))+
  #scale_colour_manual(values = col)+
  theme(panel.background = element_blank(),
        legend.key.size = unit(5, "mm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size =11,colour = "black", angle=90,hjust=.5), 
        axis.text.x = element_text(size =11,colour = "black"), 
        axis.title = element_text(size = 16, colour = "black", face="bold"),
        legend.key=element_blank(),legend.position = "right",
        legend.title = element_text(size = 12, colour = "black",face = "bold"), 
        legend.text = element_text(size = 10, colour ="black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5))+
  guides(shape = guide_legend(override.aes = list(size = 6)))+ylim(-0.5,0.5)+
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2");NMDSTrtSite

NMDSTrtSite2<-NMDSTrtSite + facet_wrap(~SampleType)+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=12, face = "bold"));NMDSTrtSite2


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
pdf("Analysis/Diversity/Graphs/NMDS-SampleName-OFT-Bac-CNF9WT2.pdf", width=8, height=6)
NMDSTrtSite
dev.off()

pdf("Analysis/Diversity/Graphs/NMDS-SampleName-OFT-Facet-Bac-CNF9WT2-legend.pdf", width=14, height=10)
NMDSTrtSite2
dev.off()


#Treatment effects & facewrap info...........................................................
pdf("Analysis/Diversity/Graphs/NMDS-Sites-TrT-OFT-Bac-CNF9WT2.pdf", width=8, height=6)
NMDSTrt
dev.off()

pdf("Analysis/Diversity/Graphs/NMDS-TrT-facetTP-OFT-Bac-CNF9WT2.pdf", width=12, height=8)
NMDSTrt2
dev.off()




#----
#----
#****************************************************************************************************----
#- .....Calculate significance ........................................................
##***************************************************************************************************----
attach(MetaRareBac)

Adon<-adonis2(OtuRareNZ~ SampleType, method = "bray",permutations = 9999,
              data=MetaRareBac);Adon # R2= 0.01721, P= 0.9979 

#Do contrast only is adonis is significant..............................................
#AdonPW<-pairwise.adonis2(OtuRareNZ ~Treatment2, data = MetaRareFun);AdonPW
capture.output(Adon,file="Analysis/Diversity/Tables/Adonis-Thawed-Original-Bac-CNF9WT2.csv")


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
MetaO<-subset(MetaRareBac,SampleType=="Original");head(MetaO[,1:5])
MetaT<-subset(MetaRareBac,SampleType=="Thawed");head(MetaT[,1:5])
MetaF<-subset(MetaRareBac,SampleType=="Frozen");head(MetaF[,1:5])
  
#Statistical significance.........................................................................
AdonO<-adonis2(OtuRareNZO~ Treatment, method = "bray",permutations = 9999,
              data=MetaO);AdonO #R2=0.23451, P= 2e-04 

AdonT<-adonis2(OtuRareNZT~ Treatment, method = "bray",permutations = 9999,
               data=MetaT);AdonT #R2=0.24035, P=2e-04 

AdonF<-adonis2(OtuRareNZF~ Treatment, method = "bray",permutations = 9999,
               data=MetaF);AdonF#R2=0.199, P=8e-04 


