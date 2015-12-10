#This code will calculate Fst between Row Type

rm(list=ls())
library(gtools)
library(hierfstat)

##INPUT FILES: from Datasets folder
GENOTYPES_POLY <-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)
Samples_info_all<-read.table("samples_information.txt",header=T)

##FORMAT FILES:
#Select sample info for those that passed QC
Samples_info<-Samples_info_all[(Samples_info_all$Alias %in% row.names(GENOTYPES_POLY)),]
#Sort by Alias
Samples_info_or<-Samples_info[order(Samples_info$Alias),]

#Convert AB genotypes to AA=1, BB=2,AB=NA
ConvertABto12<-function(dat){
	dat[dat == 'AA'] <-"1"
	dat[dat == 'BB'] <-"2"
	dat[dat =='AB'] <-"NA"
	
	return(dat)
}

Genotypes_numeric2<-as.data.frame(apply(GENOTYPES_POLY,2, ConvertABto12))


Genotypes_numeric2[1:10,1:10]
#Get SNPs names
SNPs<-colnames(Genotypes_numeric2)

#Sort genotypes by samples names
Genotypes_numeric2_or<-Genotypes_numeric2[order(row.names(Genotypes_numeric2)),]

#Put together samples info and Genotypes 
if (identical(as.character(Samples_info_or$Alias), as.character(row.names(Genotypes_numeric2_or))) == FALSE) stop (print ("Samples are in different order"))

INPUT<-cbind(as.data.frame(Samples_info_or),as.data.frame(Genotypes_numeric2_or))
#Count how many samples are in each breeding program
INPUT[1:10,1:10]
table(INPUT$Program)

TWO<-subset(INPUT,INPUT$Row_type == '2')
SIX<-subset(INPUT,INPUT$Row_type == '6')

FOCAL<-TWO[,-c(1:5)]
TOTAL<-SIX[,-c(1:5)]

##COMPUTE FST:
#Get all genotypes together, assign the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst_two_six<- cbind(as.data.frame(SNPs),Fst)
	
write.table(Fst_two_six,"Fst_RT.txt",quote=F,row.names=F,col.names=F,sep="\t")


