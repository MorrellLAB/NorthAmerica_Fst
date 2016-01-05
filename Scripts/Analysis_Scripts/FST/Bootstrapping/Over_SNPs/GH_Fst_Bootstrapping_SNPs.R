#Title: GH_Fst_Bootstrapping_SNPs.R
#Description: Count the number of outlier Fst SNPs in each iteration and how many of those
#			  are outliers when using all SNPs.
#Author: Ana Poets
#===========================================================================


rm(list=ls())
library(gtools)
library(hierfstat)

##INPUT FILES:
GENOTYPES_POLY <-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)
Samples_info_all<-read.table("samples_information.txt",header=T)


##FORMATTING FILES:
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

#Sort genotypes by samples names
Genotypes_numeric2_or<-Genotypes_numeric2[order(row.names(Genotypes_numeric2)),]

#Put together samples info and Genotypes 
if (identical(as.character(Samples_info_or$Alias), as.character(row.names(Genotypes_numeric2_or))) == FALSE) stop (print ("Samples are in different order"))

INPUT<-cbind(as.data.frame(Samples_info_or),as.data.frame(Genotypes_numeric2_or))
#Count how many samples are in each breeding program
INPUT[1:10,1:10]
table(INPUT$Program)

SPRING<-subset(INPUT,INPUT$Growth_habit == 'S')
WINTER<-subset(INPUT,INPUT$Growth_habit == 'W')

FOCAL<-SPRING[,-c(1:5)]
TOTAL<-WINTER[,-c(1:5)]

#Get all genotypes together, assign the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	
	
	
## FST BOOTSTRAPPING leaving 20% SNPs out. 100 iterations

LEAVE_OUT<-0.20

for (b in 1:100){
	LOCI_OUT<-sample(c(1:(dim(loci)[2])),round((dim(loci)[2])*(LEAVE_OUT)),replace=F)
	
	Loci_fst<-loci[, -LOCI_OUT]
	
	#Get SNPs names
	SNPs<-colnames(Loci_fst)

	#look at the Fst at each locus
		Fst<-NULL
			for (i in 1:dim(Loci_fst)[[2]]) {
			
	         Fst[i]<-varcomp(data.frame(levels,as.numeric(Loci_fst[,i])),diploid=FALSE)$F[1,1]   
	        }
	     Fst_all_GH<- cbind(as.data.frame(SNPs),Fst)
		
	write.table(Fst_all_GH,paste("Fst_GH_",b,".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
}