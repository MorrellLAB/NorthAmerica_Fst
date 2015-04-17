#This script will remove SNPs and Samples exceeding 10% missingess or heterozygosity

# remove all existing files
rm(list = ls())

#input file rows are Marker names and colums are samples names

GENOTYPES<-read.table("~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/SNP_ALLCAPs_3525_Cleaned.txt",header=T,row.names=1)
GENOTYPES[1:10,1:10]

Samples_info<-GENOTYPES[,c(1:4)]
GENOTYPES<-as.data.frame(GENOTYPES[,-c(1:4)])

# size of data set

dim(GENOTYPES)

#Sample_names <- rownames(GENOTYPES)
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

# routine to identify heterozygous SNPs
hets <- function(dat) {
dat_size <- length(na.omit(dat))
het <- which(dat == "AB")
het <- length(het)
hets_locus <- het/dat_size
return(hets_locus)
}

# all of the functions below are run using apply
# remove SNPs with missing data ≥ 10% 
# remove SNPs with observed heterozygosity > 10%
# remove SNPs with diversity = 0%; MAF = 0
miss_rm <- as.vector(which(apply(GENOTYPES,2,missing) >= 0.10))
het_rm <- as.vector(which(apply(GENOTYPES,2,hets) >= 0.10))
mono_rm <- as.vector(which(apply(GENOTYPES,2,mono) == 0))

length(mono_rm)


QC_GENOTYPES <- GENOTYPES[,c(-miss_rm,-het_rm,-mono_rm)]

#Turn table
t_QC_GENOTYPES <-as.data.frame(t(QC_GENOTYPES))
# remove lines with missing data ≥ 10% 
# remove lines with observed heterozygosity > 10%
miss_rm_accessions <- as.vector(which(apply(t_QC_GENOTYPES,2,missing) >= 0.10))
het_rm_accessions <- as.vector(which(apply(t_QC_GENOTYPES,2,hets) >= 0.10))
length(miss_rm_accessions)
# now write NA values where there 
# het_count <- apply(alch,1,hets)

QC_snp_sample_GENOTYPES <- t_QC_GENOTYPES[,c(-miss_rm_accessions,-het_rm_accessions)]

dim(QC_snp_sample_GENOTYPES)

#turn the table back to original format col.names=samples, row.names=markers

final_QC_NORTHam <-as.data.frame(t(QC_snp_sample_GENOTYPES))

#Add samples information
head(Samples_info)
samples_info_QC<-Samples_info[(row.names(Samples_info) %in% row.names(final_QC_NORTHam)),]
dim(samples_info_QC)

NORTH_AMERICAN_QC_samples_snp<-cbind(as.data.frame(samples_info_QC),as.data.frame(final_QC_NORTHam))

dim(NORTH_AMERICAN_QC_samples_snp)
#Write new table into a file
write.table(NORTH_AMERICAN_QC_samples_snp, file="~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/Data_filter_missing_hete/NORTH_AMERICAN_QC_3453samples_2396snp.txt",quote=F,sep="\t",row.names=T,col.names=T)
