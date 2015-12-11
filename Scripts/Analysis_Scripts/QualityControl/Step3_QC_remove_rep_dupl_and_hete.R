#Title           :Step3_QC_remove_rep_dupl_and_hete.R 
#Description     :Remove duplicated genotypes (<1% difference) and accessions with 
#				  >6.25% heterozygosity.
#Author		 	 :A. Poets 
#Note		     :Input file has new aliases names
#========================================================================================


# Select samples with more genotypes out of each group of duplicates
# Print tables showing the samples that needed to be removed due to duplication
# Remove accessions with >6.25% heterozygosity

rm(list=ls())
library(reshape2)
#Call table with samples count genotypes.

##INPUT FILE: Input files are in Datasets directory unless stated otherwise

#Input file is the genotypes with new aliases names after running QC Step2
NorthAm<-read.table("Barley_NorthAm_QC_AB.txt",header=T,row.names=1)

#Other samples to remove due to missing row-type or growth habit information, and the only VT Spring 6row
OTHERlines<-read.table("Lines_missing_sampleInfo.txt")
#Call the table that contain samples that have 99% similarity. So, this is a subset of the samples in that breeding program
POP<-read.csv(paste("List_NIL_duplicates.txt",sep=""),sep="\t",header=F, na.strings=c("","NA"))

Samples_info<-read.table("sample_information.txt",sep=",",header=T)

GenMap <- read.table("GeneticMap_T3_020315",header=T)


#Count how many AA and BB are in each sample. AB are treated as missing data
Count_genotypes<-function(dat) {
	AA<-length(which(dat == 'AA'))
	BB<-length(which(dat == 'BB'))
	total<-AA+BB
	}

genotypes_count<-as.data.frame(apply(NorthAm,1, Count_genotypes))

write.table(genotypes_count,"Count_genotypeCalls.txt",quote=F,row.names=T,col.names=F,sep="\t")

samples_count<-read.table("Count_genotypeCalls.txt")

#for each item in a row, find the duplicate accessions and choose the one with more data.

IdentifyUnique<-function(dat){
	position_missing<-which(is.na(dat))
	tempSamp<-rep(NA,(length(dat)-(length(position_missing))))

	for (i in 1:(length(dat)-(length(position_missing)))) {
		tempSamp[i]<-which(samples_count[,1]==dat[i])
	}
	samples_count_comparison<-samples_count[tempSamp,]
	#Select name of sample with more data 
	MoreData<-which(samples_count_comparison [,2] == max(samples_count_comparison[,2]))
	#Just in case there are two or more samples with same amount of genotypes, keep only the first entry
	To_keep<-samples_count_comparison[MoreData[1],1]
	
	}
	
Samples_to_keep<-apply(POP,1, IdentifyUnique)
Samples_to_keep_uniq<-unique(sort(Samples_to_keep))
#make a list of all the duplicated samples, so we can select those that have to be removed
All_SAMPLES<-melt(POP,id.var=colnames(POP[,1]))

Samples_to_remove<-setdiff(All_SAMPLES[,2], Samples_to_keep_uniq)
Samples_to_remove_unique<-unique(sort(Samples_to_remove))


List_samples_to_remove <-append(Samples_to_remove_unique, t(OTHERlines))
require(reshape2)
Samples_total<-length(unique(sort(All_SAMPLES[,2])))
Samples_keep_number<-length(Samples_to_keep_uniq)


#Make a table for the samples that will be removed
All_samples_to_remove_unique<-unique(List_samples_to_remove)


#Remove duplicated samples from Genotypes

NortAm_noDuplic<-subset(NorthAm,!(row.names(NorthAm) %in% All_samples_to_remove_unique))

##====Remove samples with more than 6.25% Heterozygocity (which is expected for F4 that was the minimum of inbreeding in the CAP data)

HETE<-function(dat){
	AB<-length(which(dat == 'AB'))
	TOTAL<-length(dat)
	hete_percentaje<-AB/TOTAL
}

PERCENTAGE_HETE<-apply(NortAm_noDuplic,1,HETE)

Samples_hete<-which(PERCENTAGE_HETE >0.0625)

NortAm_noDuplic_noHete<-NortAm_noDuplic[-c(Samples_hete),]

Hete_to_remove<-row.names(NortAm_noDuplic[(Samples_hete),])

All_samples_to_remove_unique_andHete<-c(All_samples_to_remove_unique, Hete_to_remove)

Sample_info_removed<-Samples_info[(Samples_info[,1] %in% All_samples_to_remove_unique_andHete),]

#List of all the samples removed from dataset
write.table(Sample_info_removed,"Table_samples_removed_duplication.txt",quote=F,row.names=F,col.names=F,sep="\t")

##### Order SNPs by Genetic Map position

DATA<-NortAm_noDuplic_noHete

if(FALSE){
	#Select lines and SNPs in AB for ACTG
	ACTG<-read.table("Barley_NorthAm_QC_ACTG.txt",header=T,row.names=1)
	Colnames<-ACTG[,(colnames(ACTG) %in% colnames(DATA))]
	Rownames<-Colnames[(row.names(Colnames) %in% row.names(DATA)),]
	dim(Rownames)
	write.table(Rownames,"Barley_NorthAm_QC_ACTG_no_duplicates.txt",quote=F,row.names=T,col.names=T,sep="\t")
	}

# 1. Select SNPs that are in our QC dataset. 2. Order by SNP name in both files 3. Turn DATA matrix & order SNPs by name 4. Add SNP genmap information 5. Order SNPs by index. 6. Remove extra columns and turn DATA matrix back again to have SNPs in columns.

SNPs_used<-intersect(colnames(DATA),GenMap$SNP)


# 1) & 2)
GenMap_used<-GenMap[(GenMap$SNP %in% SNPs_used),]
GenMap_used_or<-GenMap_used[order(GenMap_used$SNP),]

# 3)
t_DATA<-as.data.frame(t(DATA))
t_DATA_or<-t_DATA[order(row.names(t_DATA)),]

# 4)

t_DATA_map<-cbind(as.data.frame(GenMap_used_or),as.data.frame(t_DATA_or))
row.names(t_DATA_map) <-row.names(t_DATA_or)

#Check point: verify SNPs are in the same order
if (identical(row.names(t_DATA_map), as.character(t_DATA_map[,1] )) == FALSE) stop (print ("Error! SNPs are in different order"))

# 5)

t_DATA_map_or<-t_DATA_map[order(t_DATA_map$Index),]

# 6)
t_DATA_map_genotypes<-t_DATA_map_or[,-c(1:5)]

DATA_right<-as.data.frame(t(t_DATA_map_genotypes))
#Genotyping data after QC removing monomorphic, >25% missing data in  markers and accessions and >6.25% heterozygote accessions.
write.table(DATA_right,"Barley_NorthAm_QC_AB_no_duplicates.txt",quote=F,row.names=T,col.names=T,sep="\t")
