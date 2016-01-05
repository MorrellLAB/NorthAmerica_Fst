#Title: Fst_winter_spring_bootstrap.R
#Description: Fst bootstrapping- Calculate Fst between Growth habits
#            leaving out 20% of samples in each partition. 100 times.
#Author: Ana Poets
#===========================================================================

rm(list=ls())
library(gtools)
library(hierfstat)

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
SPRING<-subset(INPUT,INPUT$Growth_habit == 'S')
WINTER<-subset(INPUT,INPUT$Growth_habit == 'W')

#Bootstrapping: Randomly select 20% of samples to leave out from each partition, 100 iterations
for (b in 1:100){

LEAVE_OUT<-0.20
FOCAL<-SPRING[sample(c(1:(dim(SPRING)[1])),round((dim(SPRING)[1])*(1-LEAVE_OUT)),replace=F),-c(1:5)]
TOTAL<-WINTER[sample(c(1:(dim(WINTER)[1])),round((dim(WINTER)[1])*(1-LEAVE_OUT)),replace=F),-c(1:5)]


#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst_winter_spring<- cbind(as.data.frame(SNPs),Fst)

	write.table(Fst_winter_spring,paste("Fst_winter_spring_",b,".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
}
