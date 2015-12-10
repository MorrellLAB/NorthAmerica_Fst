#This code will calculate Fst among breeding populations
rm(list=ls())
library(gtools)
library(hierfstat)

##INPUT FILES: found in Datasets folder

GENOTYPES_POLY <-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)
Samples_info_all<-read.table("samples_information.txt",header=T)

## FILE FORMATING:
#Select sample info for those that passed QC
Samples_info<-Samples_info_all[(Samples_info_all$Alias %in% row.names(GENOTYPES_POLY)),]
dim(Samples_info)

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
GENOTYPES_ONLY<-INPUT[,-c(1:5)]


#Select the position of each program by row-type in the INPUT file
AB2<-subset(INPUT,INPUT[,3] == 'AB' & INPUT[,4]== '2')
AB6<-subset(INPUT,INPUT[,3] == 'AB' & INPUT[,4]== '6')

BA2<-subset(INPUT,INPUT[,3] == 'BA' & INPUT[,4]== '2')
BA6<-subset(INPUT,INPUT[,3] == 'BA' & INPUT[,4]== '6')
BAI2<-subset(INPUT,INPUT[,3] == 'BAI' & INPUT[,4]== '2')

MN6<-subset(INPUT,INPUT[,3] == 'MN' & INPUT[,4]== '6')

MT2<-subset(INPUT,INPUT[,3] == 'MT' & INPUT[,4]== '2')

N2<-subset(INPUT,INPUT[,3] == 'N2' & INPUT[,4]== '2')
N6<-subset(INPUT,INPUT[,3] == 'N6' & INPUT[,4]== '6')

OR2<-subset(INPUT,INPUT[,3] == 'OR' & INPUT[,4]== '2')
OR6<-subset(INPUT,INPUT[,3] == 'OR' & INPUT[,4]== '6')


UT2<-subset(INPUT,INPUT[,3] == 'UT' & INPUT[,4]== '2')
UT6<-subset(INPUT,INPUT[,3] == 'UT' & INPUT[,4]== '6')

VT6<-subset(INPUT,INPUT[,3] == 'VT' & INPUT[,4]== '6')

WA2<-subset(INPUT,INPUT[,3] == 'WA' & INPUT[,4]== '2')
WA6<-subset(INPUT,INPUT[,3] == 'WA' & INPUT[,4]== '6')

#==Select genotypes for each breeding program, 15 total=======================
BP1<-AB2[,-c(1:5)]
BP2<-AB6[,-c(1:5)]
BP3<-BA2[,-c(1:5)]
BP4<-BA6[,-c(1:5)]
BP5<-MN6[,-c(1:5)]
BP6<-MT2[,-c(1:5)]
BP7<-N2[,-c(1:5)]
BP8<-N6[,-c(1:5)]
BP9<-OR2[,-c(1:5)]
BP10<-OR6[,-c(1:5)]
BP11<-UT2[,-c(1:5)]
BP12<-UT6[,-c(1:5)]
BP13<-VT6[,-c(1:5)]
BP14<-WA2[,-c(1:5)]
BP15<-WA6[,-c(1:5)]
BP16<-BAI2[,-c(1:5)]



#Get all genotypes together
	loci<-rbind(as.data.frame(BP1),as.data.frame(BP2),as.data.frame(BP3),as.data.frame(BP4),as.data.frame(BP5),as.data.frame(BP6),as.data.frame(BP7),as.data.frame(BP8),as.data.frame(BP9),as.data.frame(BP10),as.data.frame(BP11),as.data.frame(BP12),as.data.frame(BP13),as.data.frame(BP14),as.data.frame(BP15),as.data.frame(BP16))

	levels<-c(rep(0,dim(BP1)[1]),rep(1,dim(BP2)[1]),rep(2,dim(BP3)[1]),rep(3,dim(BP4)[1]),rep(4,dim(BP5)[1]),rep(5,dim(BP6)[1]),rep(6,dim(BP7)[1]),rep(7,dim(BP8)[1]),rep(8,dim(BP9)[1]),rep(9,dim(BP10)[1]),rep(10,dim(BP11)[1]),rep(11,dim(BP12)[1]),rep(12,dim(BP13)[1]),rep(13,dim(BP14)[1]),rep(14,dim(BP15)[1]),rep(15,dim(BP16)[1]))
	head(levels)
	
## FST CALCULATION BY SNP:
	#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
     Fst_all_BP<- cbind(as.data.frame(SNPs),Fst)
	 write.table(Fst_all_BP,"Fst_BP.txt",quote=F,row.names=F,col.names=F,sep="\t")


