#This code calculate Fst comparing all populations to each other in a pairwise manner

rm(list=ls())
library(gtools)
library(hierfstat)

##INPUT FILES: from Datasets folder
GENOTYPES_POLY1 <-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)

Samples_info_all<-read.table("samples_information.txt",header=T)

## FILE FORMATTING:
#Select sample info for Samples that passed the QC

Samples_info<-Samples_info_all[(Samples_info_all$Alias %in% row.names(GENOTYPES_POLY1)),]


#Convert AB genotypes to AA=1, BB=2,AB=NA

ConvertABto12<-function(dat){
	dat[dat == 'AA'] <-"1"
	dat[dat == 'BB'] <-"2"
	dat[dat =='AB'] <-NA
	
	return(dat)
	}

Genotypes_numeric2<-as.data.frame(apply(GENOTYPES_POLY1,2, ConvertABto12))

#Get SNPs names
SNPs<-colnames(Genotypes_numeric2)

#Put together samples info and Genotypes 
#Order genotypes and samples by sample name (Alias)
Genotype_or<-Genotypes_numeric2[order(row.names(Genotypes_numeric2)),]
Samples_info_or<-Samples_info[order(Samples_info$Alias),]
if (identical(as.character(row.names(Genotype_or)),as.character(Samples_info_or[,1])) == FALSE) stop (print ("Samples are in different order"))

INPUT<-cbind(as.data.frame(Samples_info_or),as.data.frame(Genotype_or))
#Count how many samples are in each breeding program
INPUT[1:10,1:10]
table(INPUT$Program)

GENOTYPES_POLY<-INPUT[,-c(1:5)]
row.names(GENOTYPES_POLY)<-INPUT$Alias


#Select the position of each program by row-type in the INPUT file
AB2<-which(INPUT[,3] == 'AB' & INPUT[,4]== '2')
AB6<-which(INPUT[,3] == 'AB' & INPUT[,4]== '6')

BA2<-which(INPUT[,3] == 'BA' & INPUT[,4]== '2')
BA6<-which(INPUT[,3] == 'BA' & INPUT[,4]== '6')
BAI2<-which(INPUT[,3] == 'BAI' & INPUT[,4]== '2')

MN6<-which(INPUT[,3] == 'MN' & INPUT[,4]== '6')

MT2<-which(INPUT[,3] == 'MT' & INPUT[,4]== '2')

N2<-which(INPUT[,3] == 'N2' & INPUT[,4]== '2')
N6<-which(INPUT[,3] == 'N6' & INPUT[,4]== '6')

OR2<-which(INPUT[,3] == 'OR' & INPUT[,4]== '2')
OR6<-which(INPUT[,3] == 'OR' & INPUT[,4]== '6')


UT2<-which(INPUT[,3] == 'UT' & INPUT[,4]== '2')
UT6<-which(INPUT[,3] == 'UT' & INPUT[,4]== '6')

VT6<-which(INPUT[,3] == 'VT' & INPUT[,4]== '6')

WA2<-which(INPUT[,3] == 'WA' & INPUT[,4]== '2')
WA6<-which(INPUT[,3] == 'WA' & INPUT[,4]== '6')


## COMPUTING PAIRWISE FST:

POPULATIONS<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

Fst_matrix<-matrix(NA,nrow=16,ncol=16)
colnames(Fst_matrix)<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")
rownames(Fst_matrix)<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")
	
for (p in 1:length(POPULATIONS)){
	focal_pop<-get(POPULATIONS[p])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	
	#identify SNPs that in this partition are all missing data.
	Focal_NA<-which((sapply(FOCAL, function(x)all(is.na(x)))) == TRUE)
	
	#Each pop to all others
	for(o in p:length(POPULATIONS)){
		other_pop<-get(POPULATIONS[o])
		TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
		#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	
		TOTAL_NA<-which((sapply(TOTAL, function(x)all(is.na(x)))) == TRUE)
	
		# Combine both datasets and remove any SNPs that were all NA in either population
		loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
	
		ALL_snp_NA<-c(Focal_NA, TOTAL_NA)
	
		# Remove loci that are all NA
	
		if(length(ALL_snp_NA) >0){
			loci_compare<-loci[,-c(ALL_snp_NA)]} else {loci_compare<-loci}
			levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
		
			#look at the Fst at each locus

     		Fst_matrix[p,o]<-varcomp.glob(levels, loci_compare,diploid=FALSE)$F[1,1]   
        }

}	

  write.table(Fst_matrix,"Fst_matrix_ALL_pairwise.txt",quote=F,row.names=T,col.names=T,sep="\t")

