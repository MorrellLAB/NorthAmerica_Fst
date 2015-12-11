#Title           :TreeMix_input_phylogeny_w_outgroup.R 
#Description     :Format input to estimate the relationship among the data set 
#				  (spring,winter, two, six rows) using treeMix  TreeMix (Pickrell and Pritchard 2012)
#Author		 	 :A. Poets 
#Note		     :Required diploid genotyping data and an outgroup
#========================================================================================


rm(list=ls())

##INPUT FILES: from Datasets directory
NorthAm_all<-read.table("Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)

#Load samples information
samples_info<-read.table("samples_information.txt",header=T)


#Select samples for which we have genotyping data after QC
samples_info_geno<-samples_info[(samples_info$Alias %in% row.names(NorthAm_all)),]


###==== Adding an outgroup

#Genotypes from landraces downloaded from: https://github.com/AnaPoets/BarleyLandraces/blob/master/Datasets/Land_6152_SNPs_AB.txt
Landraces<-read.table("~/Dropbox/anita_301213/THESIS/Data/DATASETS_Land_wbdc_031815/Diploid/Land_6152_SNPs_AB.txt",header=T,row.names=1)

#Get SNPs that are in North America Breeding programs

SNPsharedBP<-intersect(colnames(NorthAm_all),colnames(Landraces))

#For each data set (North Am and landraces) Select Markers that are shared

NorthAm<-NorthAm_all[,(colnames(NorthAm_all) %in% SNPsharedBP)]
LANDRACE <-Landraces[,(colnames(Landraces) %in% SNPsharedBP)]

###Count allele frequency at each SNP for each population
#Using diploid unphased data
#Snps are in columns and samples in rows. The output will have SNPs in ROWs.

CountPop_landraces <-matrix(NA,nrow = (dim(land_part)[2]),ncol=length(LANDRACE))
colnames(CountPop_landraces)<-c(LANDRACE)
for (L in 1:length(LANDRACE)){
	dat <- get(LANDRACE[L])
	
	for (s in 1:(dim(dat)[2])){
		AA<-length(which(dat[,s] == 'AA'))
		BB<-length(which(dat[,s] == 'BB'))
		AB<-length(which(dat[,s] == 'AB'))
		CountA<-(2*AA)+AB
		CountB<-(2*BB)+AB
		Table_count<-cbind(CountA, CountB)
		CountPop_landraces[s,L]<-paste(Table_count,collapse=",")
	}
}



#For NorthAmerican breeding programs: Separate in populations and count allele frequencies

# How many idividuals are in each partition
as.data.frame(table(samples_info_geno $Program))
SPRING<-subset(samples_info_geno, samples_info_geno $Growth_habit == 'S')
SPRING2<-subset(samples_info_geno, samples_info_geno $Growth_habit == 'S' & samples_info_geno $Row_type == '2' )
SPRING6<-subset(samples_info_geno, samples_info_geno $Growth_habit == 'S' & samples_info_geno $Row_type == '6' )

WINTER<-subset(samples_info_geno, samples_info_geno $Growth_habit == 'W')
WINTER2<-subset(samples_info_geno, samples_info_geno $Growth_habit == 'W'& samples_info_geno $Row_type == '2' )
WINTER6<-subset(samples_info_geno, samples_info_geno $Growth_habit == 'W'& samples_info_geno $Row_type == '6' )


#select samples per population according to Spring-2row, Spring-6row, Winter2-row, Winter-6row
#Spring 2-rows
AB2s<-subset(samples_info_geno, samples_info_geno$Program == 'AB' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")
BA2s<-subset(samples_info_geno, samples_info_geno$Program == 'BA' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")
BAI2s<-subset(samples_info_geno, samples_info_geno$Program == 'BAI' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")
MT2s<-subset(samples_info_geno, samples_info_geno$Program == 'MT' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")
N2s<-subset(samples_info_geno, samples_info_geno$Program == 'N2' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")
UT2s<-subset(samples_info_geno, samples_info_geno$Program == 'UT' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")
WA2s<-subset(samples_info_geno, samples_info_geno$Program == 'WA' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "S")

#Spring 6-rows
AB6s<-subset(samples_info_geno, samples_info_geno$Program == 'AB' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "S")
BA6s<-subset(samples_info_geno, samples_info_geno$Program == 'BA' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "S")
MN6s<-subset(samples_info_geno, samples_info_geno$Program == 'MN' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "S")
N6s<-subset(samples_info_geno, samples_info_geno$Program == 'N6' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "S")
UT6s<-subset(samples_info_geno, samples_info_geno$Program == 'UT' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "S")
WA6s<-subset(samples_info_geno, samples_info_geno$Program == 'WA' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "S")

#Winter 2-rows
OR2w<-subset(samples_info_geno, samples_info_geno$Program == 'OR' & samples_info_geno$Row_type == '2' & samples_info_geno$Growth_habit == "W")

#Winter 6-rows
OR6w<-subset(samples_info_geno, samples_info_geno$Program == 'OR' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "W")
VT6w<-subset(samples_info_geno, samples_info_geno$Program == 'VT' & samples_info_geno$Row_type == '6' & samples_info_geno$Growth_habit == "W")

#Get genotypes per each population

POPULATIONS<-c("AB2s", "BA2s" ,"BAI2s" ,"MT2s" ,"N2s", "UT2s", "WA2s" ,"AB6s", "BA6s", "MN6s" ,"N6s" ,"UT6s" ,"WA6s" ,"OR2w" ,"OR6w" ,"VT6w")

for (p in 1:length(POPULATIONS)){
	pop<-get(POPULATIONS[p])
	northam_pop<-NorthAm[(row.names(NorthAm) %in% pop$Alias),]
	assign(paste(POPULATIONS[p],"_NorthAm",sep=""), northam_pop)
}


##=========
###Count allele frequency at each SNP for each population
#Using diploid unphased data
#Snps are in columns and samples in rows. The output will have SNPs in ROWs.
CountPop_NorthAm <-matrix(NA,nrow = (dim(NorthAm)[2]),ncol=length(POPULATIONS))
colnames(CountPop_NorthAm)<-c(POPULATIONS)
for (p in 1:length(POPULATIONS)){
	dat <- get(paste(POPULATIONS[p],"_NorthAm",sep=""))
	
	for (s in 1:(dim(dat)[2])){
		AA<-length(which(dat[,s] == 'AA'))
		BB<-length(which(dat[,s] == 'BB'))
		AB<-length(which(dat[,s] == 'AB'))
		CountA<-(2*AA)+AB
		CountB<-(2*BB)+AB
		Table_count<-cbind(CountA, CountB)
		CountPop_NorthAm[s,p]<-paste(Table_count,collapse=",")
	}
}


#======Combine the allele count for Landrace populations with the one from North American breeding programs.WITH outgroup
Allele_count_total<-cbind(as.data.frame(CountPop_landraces),as.data.frame(CountPop_NorthAm))

colnames(Allele_count_total)<-c("Landraces"	,"Idaho(2)",	"BuschAg(2)",	"BuschAg(Int.)"	,"Montana(2)",	"NorthDakota(2)"	,"Utah(2)",	"Washington(2)"	,"Idaho(6)"	,"BuschAg(6)",	"Minnesota(6)",	"NorthDakota(6)",	"Utah(6)",	"Washington(6)"	,"Oregon(2)",	"Oregon(6)"	,"Virginia(6)")

write.table(Allele_count_total,"TreeMix_NorthAm_land_input.gz",quote=F,row.names=F,col.names=T,sep="\t")

#add SNP names
row.names(Allele_count_total_eu_m)<-colnames(NorthAm)
write.table(Allele_count_total,"TreeMix_NorthAm_land_input_SNPnames.gz",quote=F,row.names=T,col.names=T,sep="\t")

#Run in the command line. Treemix should be installed.
#$ treemix -i TreeMix_NorthAm_land_input.gz -bootstrap -k 75 -m 3 -root Landraces -o out_run
