#Title           :Plot_smartPCA_FigS2.R
#Description     :Population structure by PCA. Plot PC1 vs PC2, as shown in figure S2
#Author		 	 :A. Poets
#========================================================================================

rm(list=ls())

#Calculate percentage of variance explained

##INPUT FILE:from Datasets directory
Samples_info<-read.table("samples_information.txt",header=T)

##Files created after running SmartPCA
EVE<-read.table("NorthAm.eval")
DATA<-read.table("NorthAm.evec",header=F,row.names=1)

#Calculate proportion of variance explained by each PC
PC1_var<-round((EVE[1,]/sum(EVE))*100,digits=1)
PC2_var<-round((EVE[2,]/sum(EVE))*100,digits=1) 
PC3_var<-round((EVE[3,]/sum(EVE))*100,digits=1) 

#Select the first 5 PC
DATA<-DATA[,1:5]

##PLOT PCA
#Separate OUTPUT by 6-rows, 2-rows, spring, winter


WINTER<-subset(Samples_info, Samples_info$Growth_habit == 'W')
SPRING<-subset(Samples_info, Samples_info$Growth_habit == 'S')

#within spring and winter divide in 2 and 6 rows

W_TWO_ROW<-subset(Samples_info, Samples_info$Row_type == '2' & Samples_info$Growth_habit == 'W')
W_SIX_ROW<-subset(Samples_info, Samples_info$Row_type == '6' &Samples_info$Growth_habit == 'W')

S_TWO_ROW<-subset(Samples_info, Samples_info$Row_type == '2' & Samples_info$Growth_habit == 'S')
S_SIX_ROW<-subset(Samples_info, Samples_info$Row_type == '6' &Samples_info$Growth_habit == 'S')

PC_WINTER<-DATA[(row.names(DATA) %in% WINTER$Alias),]
PC_SPRING<-DATA[(row.names(DATA) %in% SPRING$Alias),]
PC_W_TWO_ROW<-DATA[(row.names(DATA) %in% W_TWO_ROW $Alias),]
PC_W_SIX_ROW<-DATA[(row.names(DATA) %in% W_SIX_ROW $Alias),]
PC_S_TWO_ROW<-DATA[(row.names(DATA) %in% S_TWO_ROW $Alias),]
PC_S_SIX_ROW<-DATA[(row.names(DATA) %in% S_SIX_ROW $Alias),]


#============PLOT FigureS2: PC1 vs PC2 separated by growth habit and row-type ==========
pdf("Figures_S2A.pdf",width=7,height=5)
plot(DATA[,1],DATA[,2],col="white",main="PCA North American breeding programs",xlab=paste("PC1 (", PC1_var,"%)",sep=""), ylab=paste("PC2 (", PC2_var,"%)",sep=""))
points(PC_W_TWO_ROW[,1], PC_W_TWO_ROW[,2],cex=0.6,pch=2,col='BLUE')
points(PC_W_SIX_ROW[,1], PC_W_SIX_ROW[,2],cex=0.6,pch=3,col='BLUE')

points(PC_S_TWO_ROW[,1], PC_S_TWO_ROW[,2],cex=0.6,pch=2,col='red')
points(PC_S_SIX_ROW[,1], PC_S_SIX_ROW[,2],cex=0.6,pch=3,col='RED')

legend("bottomright",col=c("white","white","black","black"),pch=c(9,9,2,3),legend=c("Spring","Winter","Two-rowed","Six-rowed"),text.col=c("Red","blue","black","black"))

dev.off()



##=======Figure S2B-Separate by breeding program==============

#Put breeding program and row-type together
Samples_info_BP<-cbind(as.data.frame(Samples_info[,1:2]),paste(Samples_info[,3], Samples_info[,4],sep=""))
head(Samples_info_BP)

AB2<-subset(Samples_info_BP, Samples_info_BP[,3] == "AB2")
AB6<-subset(Samples_info_BP, Samples_info_BP[,3] == "AB6")
BA2<-subset(Samples_info_BP, Samples_info_BP[,3] == "BA2")
BA6<-subset(Samples_info_BP, Samples_info_BP[,3] == "BA6")
BAI2<-subset(Samples_info_BP, Samples_info_BP[,3] == "BAI2")
MN6<-subset(Samples_info_BP, Samples_info_BP[,3] == "MN6")
MT2<-subset(Samples_info_BP, Samples_info_BP[,3] == "MT2")
MT6<-subset(Samples_info_BP, Samples_info_BP[,3] == "MT6")
N2<-subset(Samples_info_BP, Samples_info_BP[,3] == "N22")
N6<-subset(Samples_info_BP, Samples_info_BP[,3] == "N66")
OR2<-subset(Samples_info_BP, Samples_info_BP[,3] == "OR2")
OR6<-subset(Samples_info_BP, Samples_info_BP[,3] == "OR6")
UT2<-subset(Samples_info_BP, Samples_info_BP[,3] == "UT2")
UT6<-subset(Samples_info_BP, Samples_info_BP[,3] == "UT6")
VT6<-subset(Samples_info_BP, Samples_info_BP[,3] == "VT6")
WA2<-subset(Samples_info_BP, Samples_info_BP[,3] == "WA2")
WA6<-subset(Samples_info_BP, Samples_info_BP[,3] == "WA6")

#GET PC per individual in each population
AB2_alias<-DATA[(row.names(DATA) %in% AB2$Alias),]
AB6_alias<-DATA[(row.names(DATA) %in% AB6$Alias),]
BA2_alias<-DATA[(row.names(DATA) %in% BA2$Alias),]
BA6_alias<-DATA[(row.names(DATA) %in% BA6$Alias),]
BAI2_alias<-DATA[(row.names(DATA) %in% BAI2$Alias),]
MN6_alias<-DATA[(row.names(DATA) %in% MN6$Alias),]
MT2_alias<-DATA[(row.names(DATA) %in% MT2$Alias),]
MT6_alias<-DATA[(row.names(DATA) %in% MT6$Alias),]
N2_alias<-DATA[(row.names(DATA) %in% N2$Alias),]
N6_alias<-DATA[(row.names(DATA) %in% N6$Alias),]
OR2_alias<-DATA[(row.names(DATA) %in% OR2$Alias),]
OR6_alias<-DATA[(row.names(DATA) %in% OR6$Alias),]
UT2_alias<-DATA[(row.names(DATA) %in% UT2$Alias),]
UT6_alias<-DATA[(row.names(DATA) %in% UT6$Alias),]
VT6_alias<-DATA[(row.names(DATA) %in% VT6$Alias),]
WA2_alias<-DATA[(row.names(DATA) %in% WA2$Alias),]
WA6_alias<-DATA[(row.names(DATA) %in% WA6$Alias),]


pdf("FigureS2B.pdf",width=7,height=5)
plot(DATA[,1],DATA[,2],col="white",main="PCA North American breeding programs",xlab=paste("PC1 (", PC1_var,"%)",sep=""), ylab=paste("PC2 (", PC2_var,"%)",sep=""))
points(AB2_alias[,1], AB2_alias[,2],cex=0.6,pch=2,col='blue')
points(AB6_alias[,1], AB6_alias[,2],cex=0.6,pch=3,col='blue')
points(BA2_alias[,1], BA2_alias[,2],cex=0.6,pch=2,col='dark green')
points(BA6_alias[,1], BA6_alias[,2],cex=0.6,pch=3,col='dark green')
points(BAI2_alias[,1], BAI2_alias[,2],cex=0.6,pch=2,col='green')
points(MN6_alias[,1], MN6_alias[,2],cex=0.6,pch=3,col='maroon')
points(MT6_alias[,1], MT6_alias[,2],cex=0.6,pch=3,col='gray')
points(N2_alias[,1], N2_alias[,2],cex=0.6,pch=2,col='purple')
points(N6_alias[,1], N6_alias[,2],cex=0.6,pch=3,col='purple')
points(OR2_alias[,1], OR2_alias[,2],cex=0.6,pch=2,col='magenta')
points(OR6_alias[,1], OR6_alias[,2],cex=0.6,pch=3,col='magenta')
points(UT2_alias[,1], UT2_alias[,2],cex=0.6,pch=2,col='cyan')
points(UT6_alias[,1], UT6_alias[,2],cex=0.6,pch=3,col='cyan')
points(VT6_alias[,1], VT6_alias[,2],cex=0.6,pch=3,col='deepskyblue4')
points(WA2_alias[,1], WA2_alias[,2],cex=0.6,pch=2,col='orange')
points(WA6_alias[,1], WA6_alias[,2],cex=0.6,pch=3,col='orange')

legend("bottomright",c("Idaho","Busch Ag.","Busch Ag.(Int.)","Minnesota","Montana","North Dakota","Oregon","Utah","Virginia","Washington", "Two-row","Six-row"),col=c("blue","dark green","green","maroon","gray","purple","magenta","cyan","deepskyblue4","orange","black","black"),title=expression(bold("Breeding Programs")),cex=0.7,pch=c(rep(15,10),2,3))

dev.off()

