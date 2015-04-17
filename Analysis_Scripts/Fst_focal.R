#This code calculate Fst in different partitions of the data:
# 1. Focal Fst for 2 and 6 rows separately, and separate by spring and winter
# 2. Focal Fst one population to all others


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

####
#AC  BA  MN  MT  N2  N6  OR  UT  VT  WA 
#382 293 384 315 384 384 298 377 327 381

####

#Select the position of each program by row-type in the INPUT file
AB2<-which(INPUT[,1] == 'AB' & INPUT[,2]== '2')
AB6<-which(INPUT[,1] == 'AB' & INPUT[,2]== '6')

BA2<-which(INPUT[,1] == 'BA' & INPUT[,2]== '2')
BA6<-which(INPUT[,1] == 'BA' & INPUT[,2]== '6')

MN6<-which(INPUT[,1] == 'MN' & INPUT[,2]== '6')

MT2<-which(INPUT[,1] == 'MT' & INPUT[,2]== '2')

N2<-which(INPUT[,1] == 'N2' & INPUT[,2]== '2')
N6<-which(INPUT[,1] == 'N6' & INPUT[,2]== '6')

OR6<-which(INPUT[,1] == 'OR' & INPUT[,2]== '6')

UT2<-which(INPUT[,1] == 'UT' & INPUT[,2]== '2')
UT6<-which(INPUT[,1] == 'UT' & INPUT[,2]== '6')

VT6<-which(INPUT[,1] == 'VT' & INPUT[,2]== '6')

WA2<-which(INPUT[,1] == 'WA' & INPUT[,2]== '2')
WA6<-which(INPUT[,1] == 'WA' & INPUT[,2]== '6')



#=== 1 .focal Fst one program to all others.======== Faster if run in MSI ~/Mohsen/New_Analysis_16042015/
if (TRUE) {
	POPULATIONS<-c("AB2","AB6","BA2","BA6","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

for (p in 1:length(POPULATIONS)){

	focal_pop<-get(POPULATIONS[p])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	TOTAL<-GENOTYPES_POLY[-c(focal_pop),]

	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		#for (i in 702) {
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst2<- cbind(as.data.frame(SNPs),Fst)
assign(paste("Fst","_",POPULATIONS[p],sep=""),Fst2)
}

for (p in 1:length(POPULATIONS)){
Infile<-get(paste("Fst_",POPULATIONS[p],sep=""))
output<-paste("~/Documents/Anita/mohsen_fst/Analysis/for_Fst/Focal_one_to_all/Fst","_",POPULATIONS[p],sep="")
write.table(Infile,output,quote=F,row.names=F,col.names=F,sep="\t")
	}
}