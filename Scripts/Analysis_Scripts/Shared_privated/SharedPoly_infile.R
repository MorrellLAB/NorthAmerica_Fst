#Title           :SharedPoly_infile.R 
#Description     :Fortmat input files for sharedPoly (libsequence package) to calculate number 
# 				   of private markers at each breeding population
#Author		 	 :A. Poets 
#Note		     :Input file has new aliases names
#========================================================================================

rm(list=ls())

##INPUT FILES: from Datasets directory.
DATA<-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)
Samples_info_all<-read.table("samples_information.txt",header=T)


#Order DATA by sample name
DATA_or<-DATA[order(row.names(DATA)),]

#Select samples passed the QC
#Sort Samples info by sample ALIAS
Samples_info<-Samples_info_all[(Samples_info_all$Alias %in% row.names(DATA)),]
Samples_info_or<-Samples_info[order(Samples_info $Alias),]


#Transform Genotypes into SharedPoly format, since SNPs are biallilic we could use either the original nuclotide calls, or vary
#the genotypes transforming the B allele to C. This is AA='A' , BB='C', AB='N', NA='N'

TRANSFORM<-function(dat){
	dat[dat == 'AA'] <- 'A'
	dat[dat == 'BB'] <- 'C'
	dat[dat == 'AB'] <- 'N'
	dat[is.na(dat)] <- 'N'
	return(dat)
}

GENOTYPES<-as.data.frame(apply(DATA_or,2,TRANSFORM))

#Create two additional rows, one with ?? and another with numbers 1...nSNP

SNPnumber<-(c(1:dim(GENOTYPES)[2]))
Question<-(rep('?',dim(GENOTYPES)[2]))

HEADER_LINE<-cbind(SNPnumber,Question)

#Separate Genotypes into populations

#Combine genotype with sample information

if (identical(as.character(Samples_info_or$Alias) , as.character(row.names(GENOTYPES))) == FALSE) stop (print ("SNPs are in different order")) 

INPUT<-cbind(as.data.frame(Samples_info_or),as.data.frame(GENOTYPES))
row.names(INPUT)<-row.names(GENOTYPES)


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

#Now put each program at the top of all samples

POPULATION<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

for (p in 1:length(POPULATION)){

Group<-get(POPULATION[p])
POP<-rbind(INPUT[Group,], INPUT[-c(Group),])
#Remove additional sample info columns

POP_clean<-POP[,-c(1:5)]

SNP_SAMPLE_count<-as.data.frame(c(dim(t_pop_clean)[2],dim(t_pop_clean)[1]))

#Contruct input file for sharedPoly. Header line shows first No.samples, No.SNPs, SNP position or order, a ? for each SNP, then add the genotypeing data.
write.table(SNP_SAMPLE_count,paste("~/Desktop/",POPULATION[p],"_vs_all.txt",sep=""),quote=F,row.names=F,col.names=F,sep="")

POP_header<-(t(HEADER_LINE))
write.table(POP_header,paste("~/Desktop/",POPULATION[p],"_vs_all.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t",append=TRUE)

write.table(POP_clean,paste("~/Desktop/",POPULATION[p],"_vs_all.txt",sep=""),quote=F,row.names=T,col.names=F,sep="\t",append=TRUE)

}

#Save a list of SNP names
SNPnames<-colnames(GENOTYPES)
write.table(SNPnames,"~/Documents/Anita/mohsen_fst/Analysis/for_sharedPoly/input/SNPinorder.txt",quote=F,row.names=T,col.names=T,sep="\t")

#run sharedPoly from the compute library to obtain number of SNPs private to each population as shown in Table 1.
