#This code will calculate Fst between Spring (2 and 6 rows) and Winter barley
#1. Comparing all spring to all winter barley
# 2. Using an even data set between spring and barley (same amount of 2 and 6 rows)

rm(list=ls())
library(gtools)
library(hierfstat)

#INPUT<-read.table("~/Dropbox/Mohsen_Fst/Data/BreedingGenotypes_homo_QC.txt",header=T)
INPUT_infile<-read.table("~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/Data_filter_missing_hete/NorthAm_3453_2367snp_geneticorder.txt",header=T)
INPUT_infile[1:10,1:10]

Samples_info<-INPUT_infile[,c(1:4)]

GENOTYPES_POLY<-INPUT_infile[,-c(1:4)]

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

#Put together samples info and Genotypes 

INPUT<-cbind(as.data.frame(Samples_info),as.data.frame(Genotypes_numeric2))
#Count how many samples are in each breeding program
INPUT[1:10,1:10]
table(INPUT$BP_Code)

SPRING<-subset(INPUT,INPUT$Growth_Habit == 'spring')
WINTER<-subset(INPUT,INPUT$Growth_Habit == 'winter')

FOCAL<-SPRING[,-c(1:4)]
TOTAL<-WINTER[,-c(1:4)]

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
	write.table(Fst_winter_spring,"~/Documents/Anita/mohsen_fst/Analysis/for_Fst/Focal_spring_winter/All_winter_all_spring/Fst_winter_spring_ALL.txt",quote=F,row.names=F,col.names=F,sep="\t")




### 2. Using equivalent samples in row-type for each spring and winter
#####  Make input files for hierfstats. Separating winter/spring and 2-rows and 6-rows at the time. Change to TRUE

	
	POPULATIONS_spring_2rows<-c("AB2","BA2","MT2","N2","UT2","WA2")
	POPULATIONS_spring_6_rows<-c("AB6","BA6","MN6","N6","UT6","WA6")
	POPULATIONS_winter_6_rows<-c("OR6","VT6")
	POPULATIONS_winter2rows_spring_2rows<-c("AB2","BA2","MT2","N2","UT2","WA2")
	
#For spring 2 rows, focal Fst
if (FALSE) {
	INPUT_FILE_spring<-subset(INPUT,INPUT$Growth_Habit == 'spring')
	
	POPULATIONS_spring_2rows<-c("AB","BA","MT","N2","UT","WA")
	INPUT_spring_2_rows<-INPUT[c(AB2,BA2,MT2,N2,UT2,WA2),]
	
	for (p in 1:length(POPULATIONS_spring_2rows)){
	
	focal_pop<-subset(INPUT_spring_2_rows, INPUT_spring_2_rows$BP_Code == POPULATIONS_spring_2rows [p])
	total_pop<-subset(INPUT_spring_2_rows, INPUT_spring_2_rows$BP_Code != POPULATIONS_spring_2rows [p])
	FOCAL<-focal_pop[,-c(1:4)]
	TOTAL<-total_pop[,-c(1:4)]

	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
dim(loci)
	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
length(levels)
	#look at the Fst at each locus

	Fst<-NULL
		for (i in 1:(dim(loci)[2])) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst<-cbind(as.data.frame(SNPs),Fst)
	assign(paste("Fst","_",POPULATIONS_spring_2rows[p],sep=""),Fst)


	output<-paste("~/Documents/Anita/mohsen_fst/Analysis/for_Fst/Focal_spring_winter/spring_2row/Fst_",POPULATIONS_spring_2rows[p],"2",sep="")
	output_file<-get(paste("Fst_", POPULATIONS_spring_2rows[p],sep=""))
	write.table(output_file,output,quote=F,row.names=F,col.names=F,sep="\t")

	}
}

#For spring 6 rows
if (FALSE) {
	INPUT_FILE_spring<-subset(INPUT,INPUT$Growth_Habit == 'spring')
	
	POPULATIONS_spring_6rows<-c("AB","BA","MN","N6","UT","WA")
	INPUT_spring_6_rows<-INPUT[c(AB6,BA6,MN6,N6,UT6,WA6),]
	
	for (p in 1:length(POPULATIONS_spring_6rows)){
	focal_pop<-subset(INPUT_spring_6_rows, INPUT_spring_6_rows$BP_Code == POPULATIONS_spring_6rows [p])
	total_pop<-subset(INPUT_spring_6_rows, INPUT_spring_6_rows$BP_Code != POPULATIONS_spring_6rows [p])
	FOCAL<-focal_pop[,-c(1:4)]
	TOTAL<-total_pop[,-c(1:4)]

	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst<-cbind(as.data.frame(SNPs),Fst)
	assign(paste("Fst","_",POPULATIONS_spring_6rows[p],sep=""),Fst)


	output<-paste("~/Documents/Anita/mohsen_fst/Analysis/for_Fst/Focal_spring_winter/spring_6row/Fst_",POPULATIONS_spring_6rows[p],"6",sep="")
	output_file<-get(paste("Fst_", POPULATIONS_spring_6rows[p],sep=""))
	write.table(output_file,output,quote=F,row.names=F,col.names=F,sep="\t")

	}
}

#For Winter 6 rows
if (FALSE) {
	INPUT_FILE_winter<-subset(INPUT,INPUT$Growth_Habit == 'winter')
	
	POPULATIONS_winter_6rows<-c("OR","VT")
	INPUT_winter_6_rows<-INPUT[c(OR6,VT6),]
	
	for (p in 1:length(POPULATIONS_winter_6rows)){
	focal_pop<-subset(INPUT_winter_6_rows, INPUT_winter_6_rows$BP_Code == POPULATIONS_winter_6rows [p])
	total_pop<-subset(INPUT_winter_6_rows, INPUT_winter_6_rows$BP_Code != POPULATIONS_winter_6rows [p])
	FOCAL<-focal_pop[,-c(1:4)]
	TOTAL<-total_pop[,-c(1:4)]

	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst<-cbind(as.data.frame(SNPs),Fst)
	assign(paste("Fst","_",POPULATIONS_winter_6rows[p],sep=""),Fst)


	output<-paste("~/Documents/Anita/mohsen_fst/Analysis/for_Fst/Focal_spring_winter/winter_6row/Fst_",POPULATIONS_winter_6rows[p],"6",sep="")
	output_file<-get(paste("Fst_", POPULATIONS_winter_6rows[p],sep=""))
	write.table(output_file,output,quote=F,row.names=F,col.names=F,sep="\t")

	}
}
