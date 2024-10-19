rm(list=ls())

setwd("/Users/fabipc/Dropbox/2-UC-Riverside/1-Dissertation/1-Research/1-HolyFire/4-Freeze-thaw-Analysis/1-Lib7Thawed-NoCNF9WT2/")

#install.packages("plyr")
library(plyr )
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(vegan)


MetaRareFun<-read.csv("Fungi/Metadata/Correlations/MetaSppRareCor-All.csv")
MetaRareBac<-read.csv("Bacteria/Metadata/Correlations/Bacteria-MergedAll-Cor.csv")


#-----
#-----
#***************************************************************************************************************----------
#---------- FUNGAL FIGURES .....................................................................................
#***************************************************************************************************************----------

    
col<-c("CNF1WT2"="#ede85f","CNF1WT5"="#d3e6f0","CNF2ET2"="#80c5e8","CNF2NT5"="#32c7cf","CNF2NT8"="#113e99",
       "CNF3WT4"="#cde8ca","CNF4ST8"="#8ba688","CNF4WT3"="#308247","CNF5ST3"="#ebe2be",
       "CNF5WT4"="#d6c467","CNF7ET2"="#856c3c","CNF7ET3"="#850001","CNF7ET4"="#d41111",
       "CNF7NT8"="#c97171","CNF7ST5"="#e3c5c5","CNF8NT3"="#c2bcbc","CNF8ST4"="#857f7f", 
       "CNF9NT5"="#d9c1e6","CNF9ST8"="#8400d6","CNF9WT2"="#c97d3a")
   
attach(MetaRareFun)
FunRegOT<-ggplot(MetaRareFun,aes(x=RichnessOriginal, y=RichnessThawed),col=SampleName)+
  geom_point(size=4, aes(col=SampleName))+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = "right")+ 
  scale_colour_manual(values = col)+ ylim(0,600)+
  labs(x = "Original Species Richness", y="")+
  annotate("text", x=350,  y=0.5,label ="Pearsons = 0.89; p < 0.0001", size=4)+
  labs(title = "A", size= 45);FunRegOT


FunRegOF<-ggplot(MetaRareFun,aes(x=RichnessOriginal, y=RichnessFrozen),col=SampleName)+
  geom_point(size=4, aes(col=SampleName))+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =14, angle=90, hjust=0.9, colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = "right")+ 
  scale_colour_manual(values = col)+ ylim(0,300)+
  labs(x = "Original Species Richness",y="Frozen Species Richness")+
  annotate("text", x=350,  y=0.5,label ="Pearsons = 0.55; p = 0.01", size=4)+
  labs(tag = "B", size= 45);FunRegOF


FunRegFT<-ggplot(MetaRareFun,aes(x=RichnessThawed, y=RichnessFrozen),col=SampleName)+
  geom_point(size=4, aes(col=SampleName))+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14, angle=90, hjust=0.9, colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = "right")+ 
  scale_colour_manual(values = col)+ ylim(0,300)+
  labs(x = "Thawed Species Richness",y="Frozen Species Richness")+
  annotate("text", x=350,  y=0.5,label ="Pearsons = 0.55; p = 0.01", size=4)+
  labs(tag = "C", size= 45);FunRegFT








#-----
#-----
#***************************************************************************************************************----------
#---------- BACTERIA FIGURES .....................................................................................
#***************************************************************************************************************----------
attach(MetaRareBac)
BacRegOT<-ggplot(MetaRareBac,aes(x=RichnessOriginal, y=RichnessThawed))+
  geom_point(size=4, aes(col=SampleID),show.legend=F)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        legend.title = element_blank(),
        legend.text = element_blank(),
      legend.position = "right")+ 
  scale_colour_manual(values = col)+ 
  labs(x = "Original Species Richness",y="Thawed Species Richness")+
  annotate("text", x=750,  y=0.5,label ="Pearsons = 0.62; p = 0.005", size=4)+
  labs(title = "A", size= 45);BacRegOT

#Automatically removes the missing sample
BacRegOF<-ggplot(MetaRareBac,aes(x=RichnessOriginal, y=RichnessFrozen),col=SampleName)+
  geom_point(size=4, aes(col=SampleName))+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14, angle=90, hjust=0.9, colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = "right")+ 
  scale_colour_manual(values = col)+ ylim(0,600)+
  labs(x = "Original Species Richness",y="Frozen Species Richness")+
  annotate("text", x=750, y=0.5,label ="Pearsons = 0.63; p = 0.004", size=4)+
  labs(tag = "B", size= 45);BacRegOF


BacRegFT<-ggplot(MetaRareBac,aes(x=RichnessThawed, y=RichnessFrozen),col=SampleName)+
  geom_point(size=4, aes(col=SampleName))+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14, angle=90, hjust=0.9, colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = "right")+ 
  scale_colour_manual(values = col)+ ylim(0,600)+
  labs(x = "Thawed Species Richness",y="Frozen Species Richness")+
  annotate("text", x=600, y=0.5,label ="Pearsons = 0.70; p = 0.001", size=4)+
  labs(tag = "C", size= 45);BacRegFT



#Export the graphs..........................................................................................----
pdf("Fungi/Analysis/Diversity/Graphs/Alpha-Cor-OrigThaw-Fun.pdf", height=6, width=8,onefile=FALSE)
FunReg
dev.off()

pdf("Fungi/Analysis/Diversity/Graphs/Alpha-Cor-OrigFroz-Fun.pdf", height=6, width=8,onefile=FALSE)
FunFrozReg
dev.off()

pdf("Fungi/Analysis/Diversity/Graphs/Alpha-Cor-ThjawFroz-Fun.pdf", height=6, width=8,onefile=FALSE)
FunFT
dev.off()


pdf("Bac-Fun-AlphaCor.pdf", height=6, width=12,onefile=FALSE)
BacFunPlots
dev.off()

pdf("Bacteria/Analysis/Diversity/Graphs/Alpha-Cor-OrigThaw-Bac.pdf", height=6, width=8,onefile=FALSE)
BacReg
dev.off()

pdf("Bacteria/Analysis/Diversity/Graphs/Alpha-Cor-OrigFroz-Bac.pdf", height=6, width=8,onefile=FALSE)
BacFrozReg
dev.off()

pdf("Bacteria/Analysis/Diversity/Graphs/Alpha-Cor-FrozThaw-Bac.pdf", height=6, width=8,onefile=FALSE)
BacFTReg
dev.off()

#-----


library(patchwork)

BacAlpha<-(BacRegOT + BacRegOF + BacRegFT +plot_layout(guides = 'collect')) +
  plot_layout(ncol = 3);BacAlpha


FunAlpha<-(FunRegOT + FunRegOF + FunRegFT +plot_layout(guides = 'collect')) +
  plot_layout(ncol = 3);FunAlpha


pdf("Bacteria/Analysis/Diversity/Graphs/Alpha-Cor-Bac-panels.pdf", height=6, width=14,onefile=FALSE)
BacAlpha
dev.off()

pdf("Fungi//Analysis/Diversity/Graphs/Alpha-Cor-Fun-panels.pdf", height=6, width=14,onefile=FALSE)
FunAlpha
dev.off()



#----
#-----
#************************************************************************************************************************************----
#Correlation significance.....................................................................................
#************************************************************************************************************************************----

BacOTCor<-cor.test(MetaRareBac$RichnessOriginal,MetaRareBac$RichnessThawed, method="pearson");BacOTCor#cor 0.;p=0.004814
BacOFCor<-cor.test(MetaRareBac$RichnessOriginal,MetaRareBac$RichnessFrozen, method="pearson");BacOFCor#cor=0.63;p=0.003545
BacTFCor<-cor.test(MetaRareBac$RichnessThawed,MetaRareBac$RichnessFrozen, method="pearson");BacTFCor#cor=0.70;p=0.001

FunOTCor<-cor.test(MetaRareFun$RichnessOriginal,MetaRareFun$RichnessThawed, method="pearson");FunOTCor#cor=0.92;p=1.186e-07
FunOFCor<-cor.test(MetaRareFun$RichnessOriginal,MetaRareFun$RichnessFrozen, method="pearson");FunOFCor#0.55;p=0.01157
FunTFCor<-cor.test(MetaRareFun$RichnessThawed,MetaRareFun$RichnessFrozen, method="pearson");FunTFCor#0.65;p=0.002



capture.output(file=BacOTCor,"Bacteria/Analysis/Diversity/Tables/Bac-PearsonCorr-OrigThawSpp.csv")
capture.output(file=FunTFCor,"Fungi/Analysis/Diversity/Tables/Fun-PearsonCorr-ThawFrozSpp.csv")
capture.output(file=BacOFCor,"Bacteria/Analysis/Diversity/Tables/Bac-PearsonCorr-OrigFrozSpp.csv")

capture.output(file=FunOFCor,"Fungi/Analysis/Diversity/Tables/Fun-PearsonCorr-OrigFrozSpp.csv")
capture.output(file=FunOTCor,"Fungi/Analysis/Diversity/Tables/Fun-PearsonCorr-OrigThawSpp.csv")
capture.output(file=FunTFCor,"Fungi/Analysis/Diversity/Tables/Fun-PearsonCorr-ThawFrozSpp.csv")




#-----
#-----
#****************************************************************************************************************************----------
#----------Fungi: avg.dist bray-curtis dissimilarity matrices - median and square root transformed for T1 T2 T3 trasplant.............
#****************************************************************************************************************************----------
#Important notes before proceeding (open) ----------------------------------------
#The analysis requires that the names match in the ASV and metadata table, however the samples that were processed in 
#lib 2 have a different naming conventions, so we will edit to  make sure that the converntions and order match, This was all done 
#in excel for ease ability

#-----
#--------End of notes ------------------------------------------------------------


#Import Metadata for bacteria and fungi-------------------------------------------------------------
#Before imporing, I checked to make sure that the names were in correct format such as (CNF#NT#)
MetaRareF<-read.csv("Fungi/Metadata/MetaRareFun.csv")
MetaRareB<-read.csv("Bacteria/Metadata/MetaRareBac.csv")



#Subset to separate thawed vs original.............................................
attach(MetaRareF)
MetaFO<-MetaRareF[which(SampleType == "Original"), ];head(MetaFO[,1:5])
MetaFT<-MetaRareF[which(SampleType == "Thawed"), ];head(MetaFT[,1:5])
MetaFF<-MetaRareF[which(SampleType == "FrozenDNA"), ];head(MetaFF[,1:5])
detach(MetaRareF)

attach(MetaRareB)
MetaBO<-MetaRareB[which(SampleType == "Original"), ];head(MetaBO[,1:2])
MetaBT<-MetaRareB[which(SampleType == "Thawed"), ];head(MetaBT[,1:2])
MetaBF<-MetaRareB[which(SampleType == "FrozenDNA"), ];head(MetaBT[,1:2])
detach(MetaRareB)

#Export the tables for any future analysis.........................................
write.csv(MetaFO,"Fungi/Metadata/MetaRareFun-Original.csv")
write.csv(MetaFT,"Fungi/Metadata/MetaRareFun-Thawed.csv")
write.csv(MetaFF,"Fungi/Metadata/MetaRareFun-Frozen.csv")


write.csv(MetaBO,"Bacteria/Metadata/MetaBacOriginal.csv")
write.csv(MetaBT,"Bacteria/Metadata/MetaBacThawed.csv")
write.csv(MetaBF,"Bacteria/Metadata/MetaBacFrozen.csv")



#Import metadata .....................................................................................................
MetaFO<-read.csv("Fungi/Metadata/Correlations/MetaRareFun-Orig.csv", row.names=1);dim(MetaFO)#20X21
MetaFT<-read.csv("Fungi/Metadata/Correlations/MetaRareFun-Thaw.csv", row.names=1);dim(MetaFT)#20X21
MetaFF<-read.csv("Fungi/Metadata/Correlations/MetaRareFun-Frozen.csv", row.names=1);dim(MetaFF)#20X21

MetaBO<-read.csv("Bacteria/Metadata/MetaRareBac-Original.csv", row.names=1);dim(MetaBO)#19x10
MetaBT<-read.csv("Bacteria/Metadata/MetaRareBac-Thawed-9s.csv", row.names=1);dim(MetaBT)#19x10
MetaBF<-read.csv("Bacteria/Metadata/MetaRareBac-Frozen-9s.csv", row.names=1);dim(MetaBF)#19x10



#Import asv tables...................................................................................................
OtuTrans2FO<-read.csv("Fungi/Analysis/ASVtables/ASVcor/OtuTrans2-Orig.csv",row.names = 1); dim(OtuTrans2FO)#20X3336
OtuTrans2FT<-read.csv("Fungi/Analysis/ASVtables/ASVcor/OtuTrans2-Thaw.csv",row.names = 1);dim(OtuTrans2FT)#20x36
OtuTrans2FF<-read.csv("Fungi/Analysis/ASVtables/ASVcor/OtuTrans2-Frozen.csv",row.names = 1);dim(OtuTrans2FF)#20X336


OtuTrans2BO<-read.csv("Bacteria/Analysis/ASVtables/OtuTrans2Cor-BacOrig.csv",row.names = 1);dim(OtuTrans2BO)#19x9429
OtuTrans2BT<-read.csv("Bacteria/Analysis/ASVtables/OtuTrans2Cor-BacThaw-9s.csv",row.names = 1);dim(OtuTrans2BT)#19x9429
OtuTrans2BF<-read.csv("Bacteria/Analysis/ASVtables/OtuTrans2Cor-BacFroz-9s.csv",row.names = 1);dim(OtuTrans2BF)#19x9429


#----
#----
#****************************************************************************************************************************************----
#Quality contol Verify that row names match...............................................................................................
#****************************************************************************************************************************************----
#Verify that row names matchs

#Fungi...............................................
all(rownames(MetaFO)==rownames(OtuTrans2FO))
all(rownames(MetaFT)==rownames(OtuTrans2FT))
all(rownames(MetaFF)==rownames(OtuTrans2FF))

#Bacteria............................................
all(rownames(MetaBO)==rownames(OtuTrans2BO))
all(rownames(MetaBT)==rownames(OtuTrans2BT))
all(rownames(MetaBF)==rownames(OtuTrans2BF))



#----
#----
#****************************************************************************************************************************************----
#QDISTANCE MATRIX AND MANTEL TEST ..............................................................................................
#****************************************************************************************************************************************----
#FO=Fungi Original , FT= Fungi thawed
DistFunO<-avgdist(OtuTrans2FO, 5786, iterations=250, meanfun=median,transf= sqrt, dmethod="bray")
DistFunT<-avgdist(OtuTrans2FT, 5786, iterations=250, meanfun=median,transf= sqrt, dmethod="bray")  
DistFunF<-avgdist(OtuTrans2FF, 5786, iterations=250, meanfun=median,transf= sqrt, dmethod="bray")  


MantFOT<-mantel(DistFunO,DistFunT, method="pearson",permutations= 9999, na.rm=TRUE);MantFOT#MantR = 0.865; p<0.0001
MantFOF<-mantel(DistFunO,DistFunF, method="pearson",permutations= 9999, na.rm=TRUE);MantFOF#MantR = 0.73; p<0.0001
MantFTF<-mantel(DistFunT,DistFunF, method="pearson",permutations= 9999, na.rm=TRUE);MantFTF#MantR = 0.6878; p<0.0001




#Calculate bacterial distance matrix..........................................................................................
DistBacO<-avgdist(OtuTrans2BO, 6515, iterations=250, meanfun=median,transf= sqrt, dmethod="bray")
DistBacT<-avgdist(OtuTrans2BT, 6515, iterations=250, meanfun=median,transf= sqrt, dmethod="bray")
DistBacF<-avgdist(OtuTrans2BF, 6515, iterations=250, meanfun=median,transf= sqrt, dmethod="bray")

#Calculate mantel test
MantBOT<-mantel(DistBacO,DistBacT, method="pearson",permutations= 9999, na.rm=TRUE);MantBOT#MantR 0.9499; p<0.0001 
MantBOF<-mantel(DistBacO,DistBacF, method="pearson",permutations= 9999, na.rm=TRUE);MantBOF#MantR 0.8716, p<0.0001
MantBTF<-mantel(DistBacT,DistBacF, method="pearson",permutations= 9999, na.rm=TRUE);MantBTF#MantR 0.8845; p<0.0001


#----                                                                                                                                                                                                                                                      
#----
#******************************************************************************************************************************-----
#.....Convert the distance matrix into a vector so taht we can use it in ggplot......................
#******************************************************************************************************************************-----
#library(gdata)

#Fungi............................................
#Thaw......................................
FThaw<-as.data.frame(as.matrix(DistFunT));FThaw
FunThaw<-FThaw[lower.tri(FThaw)];FunThaw #introduces NA, wehre diagnol is 0


#Original..................................
FOrig<-as.data.frame(as.matrix(DistFunO));FOrig
FunOrig<-FOrig[lower.tri(FOrig)];FunOrig


#Frozen..................................
FFroz<-as.data.frame(as.matrix(DistFunF));FFroz
FunFroz<-FFroz[lower.tri(FFroz)];FunFroz



FunOTdata<- as.data.frame(cbind(FunOrig,FunThaw));FunOTdata#make it a dataframe
FunOFdata<- as.data.frame(cbind(FunOrig,FunFroz));FunOFdata#make it a dataframe
FunTFdata<- as.data.frame(cbind(FunThaw,FunFroz));FunTFdata#make it a dataframe


#Export data......................................................................................
write.csv(FunOTdata,"Fungi/Analysis/Diversity/Tables/Fungi-Dist-OriginalThaw.csv")
write.csv(FunOFdata,"Fungi/Analysis/Diversity/Tables/Fungi-Dist-OriginalFrozen.csv")
write.csv(FunTFdata,"Fungi/Analysis/Diversity/Tables/Fungi-Dist-ThawFrozen.csv")

#Write full datamatrix to make sure that the bing worked..................................
write.csv(FThaw,"Fungi/Analysis/Diversity/Tables/Fungi-Dist-Thaw.csv")
write.csv(FOrig,"Fungi/Analysis/Diversity/Tables/Fungi-Dist-Original.csv")
write.csv(FFroz,"Fungi/Analysis/Diversity/Tables/Fungi-Dist-Frozen.csv")



#Bacteria............................................................................................

#Thaw..............................................
BThaw<-as.data.frame(as.matrix(DistBacT));BThaw
BacThaw<-BThaw[lower.tri(BThaw)];BacThaw 


#Original..........................................
BOrig<-as.data.frame(as.matrix(DistBacO));BOrig
BacOrig<-BOrig[lower.tri(BOrig)];BacOrig


#Frozen............................................
BFroz<-as.data.frame(as.matrix(DistBacF));BFroz
BacFroz<-BFroz[lower.tri(BFroz)];BacFroz


#Create a list of the bind data ...................................................
BacOTdata<- as.data.frame(cbind(BacOrig,BacThaw));head(BacOTdata[1:3,])#make it a dataframe
BacOFdata<- as.data.frame(cbind(BacOrig,BacFroz));head(BacOFdata[1:3,])#make it a dataframe
BacTFdata<- as.data.frame(cbind(BacThaw,BacFroz));head(BacTFdata[1:3,])#make it a dataframe



#Write matrix..........................................................................................
write.csv(BacOTdata,"Bacteria/Analysis/Diversity/Tables/Bac-Dist-OriginalThaw.csv")
write.csv(BacOFdata,"Bacteria/Analysis/Diversity/Tables/Bac-Dist-OriginalFrozen.csv")
write.csv(BacTFdata,"Bacteria/Analysis/Diversity/Tables/Bac-Dist-ThawFrozen.csv")

#Write full datamatrix to make sure that the bing worked..................................
write.csv(BThaw,"Bacteria/Analysis/Diversity/Tables/Bac-Dist-Thaw.csv")
write.csv(BOrig,"Bacteria/Analysis/Diversity/Tables/Bac-Dist-Original.csv")
write.csv(BFroz,"Bacteria/Analysis/Diversity/Tables/Bac-Dist-Frozen.csv")




#----
#----
#********************************************************************************************************************************-----
#.....CREATE CORRELATION 1:1 GRAPHS.....................
#******************************************************************************************************************************-----
FunOTdata<-read.csv("Fungi/Analysis/Diversity/Tables/Fungi-Dist-OriginalThaw.csv", row.names = 1)
FunOFdata<-read.csv("Fungi/Analysis/Diversity/Tables/Fungi-Dist-OriginalFrozen.csv", row.names = 1)
FunTFdata<-read.csv("Fungi/Analysis/Diversity/Tables/Fungi-Dist-ThawFrozen.csv", row.names = 1)

BacOTdata<-read.csv("Bacteria/Analysis/Diversity/Tables/Bac-Dist-OriginalThaw.csv", row.names = 1)
BacOFdata<-read.csv("Bacteria/Analysis/Diversity/Tables/Bac-Dist-OriginalFrozen.csv", row.names = 1)
BacTFdata<-read.csv("Bacteria/Analysis/Diversity/Tables/Bac-Dist-ThawFrozen.csv", row.names = 1)

#Rename the columns of the datasheet .........................................................................
names(FunOTdata) <- c("SampleID","Original","Thawed");head(FunOTdata[1:3])
names(FunOFdata) <- c("SampleID","Original","Frozen");head(FunOFdata[1:3])
names(FunTFdata) <- c("SampleID","Thawed","Frozen");head(FunTFdata[1:3])

names(BacOTdata) <- c("SampleID","Original","Thawed");head(BacOTdata[1:3])
names(BacOFdata) <- c("SampleID","Original","Frozen");head(BacOFdata[1:3])
names(BacTFdata) <- c("SampleID","Thawed","Frozen");head(BacTFdata[1:3])


#-----
#---------
#CREATE FUNGAL GRAPHS...........................................................................................----
col<-c("CNF1WT5"="#d3e6f0","CNF2ET2"="#80c5e8","CNF2NT5"="#32c7cf","CNF2NT8"="#113e99",
       "CNF3WT4"="#cde8ca","CNF4ST8"="#8ba688","CNF4WT3"="#308247","CNF5ST3"="#ebe2be",
       "CNF5WT4"="#d6c467","CNF7ET2"="#856c3c","CNF7ET3"="#850001","CNF7ET4"="#d41111",
       "CNF7NT8"="#c97171","CNF7ST5"="#e3c5c5","CNF8NT3"="#c2bcbc","CNF8ST4"="#857f7f", 
       "CNF9NT5"="#d9c1e6","CNF9ST8"="#8400d6","CNF9WT2"="#c97d3a")


CorFOT<-ggplot(FunOTdata,aes(x=Thawed, y=Original))+
  geom_point(size=2.5, aes(col=SampleID),show.legend=F)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  scale_colour_manual(values = col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =13, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 13,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "right")+ 
  xlim(0.4,1.0)+labs(x = "Thaw Beta Diversity", y="Original Beta Diversity")+
  annotate("text", x=0.8,  y=0.4,label ="Mantel R = 0.87; p < 0.0001", size=3.5)+
  labs(title = "A", size= 45);CorFOT

CorFOF<-ggplot(FunOFdata,aes(x=Frozen, y=Original))+
  geom_point(size=2.5, aes(col=SampleID),show.legend=T)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  scale_colour_manual(values = col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =13, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 13,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "right")+ 
  xlim(0.4,1.0)+labs(x = "Frozen Beta Diversity", y="Original Beta Diversity")+
  annotate("text", x=0.8,  y=0.4,label ="Mantel R = 0.73; p < 0.0001", size=3.5)+
  labs(title = "B", size= 45);CorFOF


CorFTF<-ggplot(FunTFdata,aes(x=Thawed, y=Frozen))+
  geom_point(size=2.5, aes(col=SampleID),show.legend=T)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  scale_colour_manual(values = col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =13, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 13,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "right")+ 
  xlim(0.4,1.0)+labs(x = "Thaw Beta Diversity", y="Frozen Beta Diversity")+
  annotate("text", x=0.75,  y=0.4,label ="Mantel R = 0.69; p < 0.0001", size=4.5)+
  labs(title = "C", size= 45);CorFTF



#CREATE BACTERIAL GRAPHS...........................................................................................----
CorBOT<-ggplot(BacOTdata,aes(x=Thawed, y=Original))+
  geom_point(size=2.5, aes(col=SampleID),show.legend=F)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  scale_colour_manual(values = col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =13, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 13,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "right")+ 
 labs(x = "Thawed Beta Diversity", y="Original Beta Diversity")+
  annotate("text", x=0.8,  y=0.3,label ="Mantel R = 0.95; p < 0.0001", size=3.5)+
  labs(title = "A", size= 45);CorBOT

CorBOF<-ggplot(BacOFdata,aes(x=Frozen ,y=Original))+
  geom_point(size=2.5, aes(col=SampleID),show.legend=T)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  scale_colour_manual(values = col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =13, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 13,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_blank(),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "right")+ 
  labs(x = "Frozen Beta Diversity", y="Original Beta Diversity")+
  annotate("text", x=0.8,  y=0.3,label ="Mantel R = 0.87; p < 0.0001", size=3.5)+
  labs(title = "B", size= 45);CorBOF


CorBTF<-ggplot(BacTFdata,aes(x=Thawed, y=Frozen))+
  geom_point(size=2.5, aes(col=SampleID),show.legend=T)+ 
  geom_smooth(method=lm, se=TRUE, alpha=0.2, col="black")+
  scale_colour_manual(values = col)+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =13, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 13,colour = "black"), 
        axis.title = element_text(size = 14, colour = "black"), 
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "right")+ 
  labs(x = "Thaw Beta Diversity", y="Frozen Beta Diversity")+
  annotate("text", x=0.75,  y=0.3,label ="Mantel R = 0.88; p < 0.0001", size=4)+
  labs(title = "C", size= 45);CorBTF




#----
#----
#****************************************************************************************************-----
#.....CREATE PANELS OF THE ALPGA AND BETA GRAPHS.....................
#****************************************************************************************************-----
#Create patchword
library(patchwork)

FungalPlots<-(CorFOT+ CorFOF + CorFTF + plot_layout(guides = 'collect')) +
  plot_layout(ncol = 3);FungalPlots

FungalPlots1<-(CorFOT+ CorFOF + plot_layout(guides = 'collect')) +
  plot_layout(ncol = 2);FungalPlots1



BacPlots<-(CorBOT+ CorBOF + CorBTF + plot_layout(guides = 'keep')) +
  plot_layout(ncol = 3);BacPlots

BacPlots1<-(CorBOT+ CorBOF + plot_layout(guides = 'keep')) +
  plot_layout(ncol = 2);BacPlots1


#Fungal plots.................................................................................................
pdf("Fungi/Analysis/Diversity/Graphs/Beta-Reg-All-Fun.pdf", height=6, width=14,onefile=FALSE)
FungalPlots
dev.off()

pdf("Fungi/Analysis/Diversity/Graphs/Beta-Reg-3-Fun.pdf", height=6, width=10,onefile=FALSE)
FungalPlots1
dev.off()


#Bacterial plots.................................................................................................
pdf("Bacteria/Analysis/Diversity/Graphs/Beta-Reg-All-Bac.pdf", height=6, width=10,onefile=FALSE)
BacPlots
dev.off()

pdf("Bacteria/Analysis/Diversity/Graphs/Beta-Reg-3-Bac.pdf", height=6, width=10,onefile=FALSE)
BacPlots1
dev.off()

