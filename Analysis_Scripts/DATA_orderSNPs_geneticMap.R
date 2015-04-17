#Find SNPs that have genetic map position. Only mapped SNPs will be used for Fst and PHS analysis

INPUT_file<-read.table("~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/Data_filter_missing_hete/NORTH_AMERICAN_QC_3453samples_2396snp.txt",header=T)
INPUT_file[1:10,1:10]
Samples_identifier<-INPUT_file[,c(1:4)]

GENOTYPES_POLY<-INPUT_file[,-c(1:4)]
GENOTYPES_POLY[1:10,1:10]
#Get SNPs names
SNPs<-colnames(GENOTYPES_POLY)

genmap<-read.table("~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/GeneticMap_T3_020315",header=T)
genmap[1:10,]

#Order genmap by SNP name
genmap_or<-genmap[order(genmap$SNP),]

#Order GENOTYPES by SNP name

GENOTYPES_or<-GENOTYPES_POLY[,order(colnames(GENOTYPES_POLY))]
GENOTYPES_or[1:10,1:10]

#Find genmap SNPs present in GENOTYPES_or

SNPs_used<-intersect(genmap_or$SNP, colnames(GENOTYPES_or))

#Select SNPs used in both data sets

GENOTYPES_snp<-GENOTYPES_or[,(colnames(GENOTYPES_or) %in% SNPs_used)]
genmap_or_snp<-genmap_or[(genmap_or$SNP %in% SNPs_used),]

#Put SNPs info and genotypes together, order SNPs by index 

t_genotypes<-as.data.frame(t(GENOTYPES_snp))
t_genotypes[1:10,1:10]

GENOTYPES_SNP_info<-cbind(as.data.frame(genmap_or_snp),as.data.frame(t_genotypes))

GENOTYPES_SNP_info_index<-GENOTYPES_SNP_info[order(GENOTYPES_SNP_info$Index),]
GENOTYPES_SNP_info_index[1:10,1:10]

justGenotypes<-GENOTYPES_SNP_info_index[,-c(2:5)]
t_justGenotypes<-as.data.frame(t(justGenotypes))
t_justGenotypes[1:10,1:10]

#Write and Read again to accomodate header
write.table(t_justGenotypes,"~/Desktop/temp.txt",quote=F,row.names=F,col.names=F,sep="\t")
DATA<-read.table("~/Desktop/temp.txt",header=T)

DATA[1:10,1:10]

#Add samples information
SNP_ALLCAPs_2396_GENMAP<-cbind(as.data.frame(Samples_identifier),as.data.frame(DATA))

write.table(SNP_ALLCAPs_2396_GENMAP,"~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/Data_filter_missing_hete/NorthAm_3453_2367snp.txt",quote=F,row.names=F,col.names=T,sep="\t")