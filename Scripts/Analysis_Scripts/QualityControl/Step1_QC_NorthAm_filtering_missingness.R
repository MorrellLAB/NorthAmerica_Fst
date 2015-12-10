#This script will remove SNPs and Samples exceeding 25% missingess and monomorphic SNPs. (Didn't use the heterozygosity filter)

# remove all existing files
rm(list = ls())

##Input file is the genotyping file downloaded from www.triticeaetoolbox.org
ORIGINAL<-read.csv("snpfile.txt",sep="\t",header=T,row.names=1,na.strings="")

Replace_dashforNA<-function(dat){
	dat[dat =='-'] <-"NA"
	return(dat)
}

GENOTYPES<-as.data.frame(apply(ORIGINAL, 2, Replace_dashforNA))

write.table(GENOTYPES,"BARLEYgenotypes_AB.txt",quote=F,row.names=T,col.names=T,sep="\t")

#input file rows are Marker's names, and colums are sample names

GENOTYPES<-read.table("BARLEYgenotypes_AB.txt",header=T,row.names=1)
SNPs<-colnames(GENOTYPES)


# routine to identify SNPs with missing data, shown as '--'
missing <- function(dat) { 
#dat <- na.omit(dat)
miss <- which(is.na(dat))
miss <- length(miss)
dat_size <- length(dat)
miss_locus <- miss/dat_size
return(miss_locus)
}

# routine to identify monomorphic SNPs
# actually identifies minor allele frequency (MAF)
mono <- function(dat) {
#dat_size <- length(na.omit(dat))
AA <- dat[dat == "AA"]
AA <- length(AA)
BB <- dat[dat == "BB"]
BB <- length(BB)
smaller <- min(c(AA,BB))
return(smaller)
}

# all of the functions below are run using apply
# remove SNPs with missing data ≥ 25% 
# remove SNPs with diversity = 0%; MAF = 0
miss_rm <- as.vector(which(apply(GENOTYPES,2,missing) >= 0.25))
mono_rm <- as.vector(which(apply(GENOTYPES,2,mono) == 0))
QC_GENOTYPES <- GENOTYPES[,c(-miss_rm,-mono_rm)]

#Turn table
t_QC_GENOTYPES <-as.data.frame(t(QC_GENOTYPES))
# remove ACCESSIONS with missing data ≥ 25% 
miss_rm_accessions <- as.vector(which(apply(t_QC_GENOTYPES,2,missing) >= 0.25))
final_QC_NORTHam <- t_QC_GENOTYPES[,c(-miss_rm_accessions)]

#Write new table into a file
write.table(final_QC_NORTHam, file="Barley_genotypes_QC_AB.txt",quote=F,sep="\t",row.names=T,col.names=T)

####==============Get the same set of samples and SNPs from the ACTG file

#Input file downloaded from www.triticeaetoolbox.org with parameters as reported in the manuscript.

BARLEY_ACTG<-read.table("genotype.hmp.txt",header=F,row.names=1)

#Turn table to have SNPs in columns and samples in rows. Export table and import again. THis will add an X infront of the SNPs, to match with the names in the QC_AB file.

BARLEY_ACTG_GEN<-BARLEY_ACTG[,-c(1:3)]
t_BARLEY_ACTG_GEN <-as.data.frame(t(BARLEY_ACTG_GEN))

write.table(t_BARLEY_ACTG_GEN,"BARLEYgenotypes_ACTG.txt",quote=F,row.names=F,col.names=T,sep="\t")

BARLEY_actg<-read.table("BARLEYgenotypes_ACTG.txt",header=T,row.names=1)

#Select samples that pass the QC

ACTG_LINES<-BARLEY_actg[(row.names(BARLEY_actg) %in% row.names(final_QC_NORTHam)),]

ACTG_SNPS<-ACTG_LINES[,(colnames(ACTG_LINES) %in% colnames(final_QC_NORTHam))]

write.table(ACTG_SNPS,"Barley_genotypes_QC_ACTG.txt",quote=F,row.names=T,col.names=T,sep="\t")
