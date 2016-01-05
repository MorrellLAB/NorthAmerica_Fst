#Title: Fst_six_two_row_bootstrap.R
#Description: Fst bootstrapping- Calculate Fst between Row Type 
#            leaving out 20% of samples in each partition. 100 times.
#Author: Ana Poets
#===========================================================================

rm(list=ls())
library(gtools)
library(hierfstat)


##INPUT FILES:
GENOTYPES_POLY <-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)
Samples_info_all<-read.table("samples_information.txt",header=T)

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

##BOOTSTRAPPING OVER SAMPLES:
#Randomly select 20% individuals from each group to leave out, 100 iterations
for (b in 1:100){

LEAVE_OUT<-0.20
FOCAL<-TWO[sample(c(1:(dim(TWO)[1])),round((dim(TWO)[1])*(1-LEAVE_OUT)),replace=F),-c(1:5)]
TOTAL<-SIX[sample(c(1:(dim(SIX)[1])),round((dim(SIX)[1])*(1-LEAVE_OUT)),replace=F),-c(1:5)]

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
	
write.table(Fst_two_six,paste("Fst_RT_",b,".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

}
