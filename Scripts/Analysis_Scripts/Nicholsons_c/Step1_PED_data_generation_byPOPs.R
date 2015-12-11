#Title           :Step1_PED_data_generation_byPOPs.R
#Description     :Convert genotype file into .PED format for Nicholson Fst input data generation. 
#			      Separating spring vs winter and 2 from 6 rows phenotypes
#Author		 	 :A. Poets
#Date			 :May 5, 2015
#Note		     :Requires diploid genotypes
#========================================================================================

rm(list=ls())

##INPUT FILES: located in Datasets directory
Sample_info_all<-read.table("samples_information.txt",header=T)
NorthAm<-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)

#Get samples that passed QC
Sample_info<-Sample_info_all[(Sample_info_all$Alias %in% row.names(NorthAm)),]

#Separate spring and winter, and 6 and 2 rows
spring6<-subset(Sample_info, Sample_info $Growth_habit == 'S' & Sample_info $Row_type == 6)
spring2<-subset(Sample_info, Sample_info $Growth_habit == 'S' & Sample_info $Row_type == 2)
winter6<-subset(Sample_info, Sample_info $Growth_habit == 'W' & Sample_info $Row_type == 6)

#Divide genotypes into the 3 groups.

INPUT_spring6<-NorthAm[(row.names(NorthAm) %in% spring6$Alias),]
INPUT_spring2<-NorthAm[(row.names(NorthAm) %in% spring2$Alias),]
INPUT_winter6<-NorthAm[(row.names(NorthAm) %in% winter6$Alias),]

#order populations by breeding program names
INPUT_spring6<-INPUT_spring6[order(row.names(INPUT_spring6)),]
INPUT_spring2 <-INPUT_spring2[order(row.names(INPUT_spring2)),]
INPUT_winter6 <-INPUT_winter6[order(row.names(INPUT_winter6)),]

PARTITIONS<-c("spring6", "spring2", "winter6")

for (i in 1:length(PARTITIONS)){
	
	GENOTYPES <-get(paste("INPUT_",PARTITIONS[i],sep=""))
	
#Removing monomorphic SNPs

MONO<-function(dat){
	AA<-length(which(dat == 'AA'))
	BB<-length(which(dat == 'BB'))
	NONE<-length(which(is.na(dat)))
	MAX<-max(AA,BB)
	if(((length(dat) - NONE)-MAX) > 0) {'poly'} else {"mono"}
	}

TOTAL_calls<-apply(GENOTYPES,2,MONO)
mono_rm<-which(TOTAL_calls == "mono")

if (length(mono_rm) >0) {
GENOTYPES<-GENOTYPES[,-c(mono_rm)] }


GENOTYPES[1:10,1:10]
#Create extra columns ID information just with 0

Extra<-rep(0,(dim(GENOTYPES)[1]))

#Formatting PED file ID columns info

ID<-cbind(Extra,as.data.frame(row.names(GENOTYPES)), Extra, Extra, Extra, Extra)
head(ID)
GENOTYPES[1:10,1:10]

#Set genotypes in the right format

genotypesFormat<-function(dat){
	dat[dat == 'AA'] <-c("A A")
	dat[dat == 'BB'] <-c("B B")
	dat[dat == 'AB'] <-c("A B")
	dat[is.na(dat)] <-c("0 0")
	return(dat)
}

Genotypes_split<-as.data.frame(t(apply(GENOTYPES, 1,genotypesFormat)))
SNPs<-colnames(Genotypes_split)

#Putting together PED file format
PED<-cbind(as.data.frame(ID),as.data.frame(Genotypes_split))

write.table(PED,paste("NorthAm_",PARTITIONS[i],".ped",sep=""),quote=F,row.names=F,col.names=F,sep=" ")
write.table(PED,paste("NorthAm_w_header_",PARTITIONS[i],".ped",sep=""),quote=F,row.names=F,col.names=T,sep=" ")
write.table(SNPs,paste("SNPs_NorthAm_nicholson_",PARTITIONS[i],".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

#Create .CLST sample ID and cluster membership


AB_2 <-subset(Sample_info, Sample_info$Program == 'AB' & Sample_info$Row_type == 2)
BA_2 <-subset(Sample_info, Sample_info$Program == 'BA' & Sample_info$Row_type == 2)
BAI_2 <-subset(Sample_info, Sample_info$Program == 'BAI' & Sample_info$Row_type == 2)
MT_2 <-subset(Sample_info, Sample_info$Program == 'MT' & Sample_info$Row_type == 2)
N_2 <-subset(Sample_info, Sample_info$Program == 'N2' & Sample_info$Row_type == 2)
OR_2 <-subset(Sample_info, Sample_info$Program == 'OR' & Sample_info$Row_type == 2)
UT_2 <-subset(Sample_info, Sample_info$Program == 'UT' & Sample_info$Row_type == 2)
WA_2 <-subset(Sample_info, Sample_info$Program == 'WA' & Sample_info$Row_type == 2)


AB_6 <-subset(Sample_info, Sample_info$Program == 'AB' & Sample_info$Row_type == 6)
BA_6 <-subset(Sample_info, Sample_info$Program == 'BA' & Sample_info$Row_type == 6)
MN_6 <-subset(Sample_info, Sample_info$Program == 'MN' & Sample_info$Row_type == 6)
N_6 <-subset(Sample_info, Sample_info$Program == 'N6' & Sample_info$Row_type == 6)
UT_6 <-subset(Sample_info, Sample_info$Program == 'UT' & Sample_info$Row_type == 6)
WA_6 <-subset(Sample_info, Sample_info$Program == 'WA' & Sample_info$Row_type == 6)
OR_6 <-subset(Sample_info, Sample_info$Program == 'OR' & Sample_info$Row_type == 6)
VT_6 <-subset(Sample_info, Sample_info$Program == 'VT' & Sample_info$Row_type == 6)


#Get the names of all the samples in each group together, along with the population classifier


if (PARTITIONS[i] == "spring6"){
NAMES_ONLY<-c(as.character(AB_6$Alias), as.character(BA_6$Alias), as.character(MN_6$Alias), as.character(N_6$Alias), as.character(UT_6$Alias), as.character(WA_6$Alias))

POP_ASSIG<-c(rep(1,dim(AB_6)[1]),rep(2,dim(BA_6)[1]),rep(3,dim(MN_6)[1]),rep(4,dim(N_6)[1]),rep(5,dim(UT_6)[1]),rep(6,dim(WA_6)[1]))}

if (PARTITIONS[i] == "spring2"){
	NAMES_ONLY<-c(as.character(AB_2$Alias),as.character(BA_2$Alias),as.character(BAI_2$Alias),as.character(MT_2$Alias),as.character(N_2$Alias),as.character(UT_2$Alias),as.character(WA_2$Alias))
	POP_ASSIG <-c(rep(1,dim(AB_2)[1]),rep(2,dim(BA_2)[1]),rep(3,dim(BAI_2)[1]),rep(4,dim(MT_2)[1]),rep(5,dim(N_2)[1]),rep(6,dim(UT_2)[1]),rep(7,dim(WA_2)[1]))
	}

if (PARTITIONS[i] == "winter6"){
	NAMES_ONLY<-c(as.character(OR_6$Alias),as.character(VT_6$Alias))
	POP_ASSIG <-c(rep(1,dim(OR_6)[1]),rep(2,dim(VT_6)[1]))
}


CLST<-cbind(as.data.frame(NAMES_ONLY),as.data.frame(POP_ASSIG))

write.table(CLST,paste("Samples_pops_",PARTITIONS[i],".clst",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
}

#Now use Step2_Make_Nicholson_FST_inputs.py to create all the input files required by Nicholson's Fst
#python Step2_Make_Nicholson_FST_Inputs.py -p NorthAm_spring2.ped -c Samples_pops_spring2.clst -o NorthAm_spring2
#Then use Step3_Nicholson_Fst_run.R 