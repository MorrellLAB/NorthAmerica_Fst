rm(list=ls())

NorthAm<-read.table("~/Documents/Anita/mohsen_fst/Data/t3_download_04272015/ANA_names_Barley_NorthAm_QC_AB_no_duplicates_or.txt",header=T,row.names=1)


#Identify how many of those SNPs have genetic position

genmap<-read.table("~/Documents/Anita/mohsen_fst/Data/t3_download_04272015/GeneticMap_T3_020315",header=T)

head(genmap)

SNPs_common_map<-intersect(colnames(NorthAm),genmap$SNP)
length(SNPs_common_map)

#Select SNPs position for those in the panel of NorthAm

genmap_northam<-genmap[(genmap$SNP %in% SNPs_common_map),]

#check that SNPs are in the same order as in the genetic map

if (identical(as.character(colnames(NorthAm)),as.character(genmap_northam $SNP)) == FALSE) stop (print ("Error! SNPs are not in genetic map order"))

#Add SNPs from genetic map and order by genetic map order using the index

t_NorthAm<-t(NorthAm)

north_map<-cbind(genmap_northam, t_NorthAm)

#Sort by Index

NorthAm_genmap_or<-north_map[order(north_map $Index),]

#Remove extra columns and turn tables to have SNPs in columnas and samples in rows

NorthAm_genmap_or_ready<-NorthAm_genmap_or[,-c(2:5)]
NorthAm_genmap_or_ready_t<-t(NorthAm_genmap_or_ready)

#Make the first row the column headers
colnames(NorthAm_genmap_or_ready_t)<-NorthAm_genmap_or_ready_t[1,]
NorthAm_genmap_or_ready_t <-as.data.frame(NorthAm_genmap_or_ready_t[-1,])
NorthAm_genmap_or_ready_t[1:10,1:10]



######################################
##Separate Genotypes AA in A A, BB in B B, AB in A B, NA in ? ?
#####Select NorthAm_genmap_or_ready_t

DATA<-NorthAm_genmap_or_ready_t

SNPnames<-colnames(DATA)
Sep_genotypes<-function(dat){
	dat [dat == 'AA'] <-c('A A')
	dat [dat == 'BB'] <-c('B B')
	dat [dat == 'AB'] <-c('A B')
	dat [is.na(dat)] <-c('? ?') 
	return (dat)
}
DOUBLED_GENOTYPES<-t(as.data.frame(apply(DATA,1,Sep_genotypes)))
write.table(DOUBLED_GENOTYPES,"~/Desktop/temp.txt",quote=F,row.names=T,col.names=F,sep="\t")

#Import file again without header so we can split the genotypes and add SNPs names

Genotypes<-read.table("~/Desktop/temp.txt",header=F,row.names=1)
SAMPLE<-row.names(Genotypes)
#Create a string of Sample names that had and extra _2, so we can merge the even and odd genotypes and get only one file 

SAMPLE_DOUBLE<-paste(SAMPLE,"_2",sep="")

#Separate genotypes in haplotypes. Add SNPs and Samples_2 names
ODD_columns<-Genotypes[,seq(1,ncol(Genotypes),by=2)]
colnames(ODD_columns) <-SNPnames
EVEN_columns<-Genotypes[,seq(2,ncol(Genotypes),by=2)]
row.names(EVEN_columns) <-SAMPLE_DOUBLE
colnames(EVEN_columns)<-SNPnames

#Combine both files

DOUBLED_DATA<-rbind(as.data.frame(ODD_columns), as.data.frame(EVEN_columns))

#Sort by sample name

DOUBLED_DATA_or <-DOUBLED_DATA[order(row.names(DOUBLED_DATA)),]

#Devide samples by chromosome
Chr1<-subset(genmap_northam, genmap_northam $chromosome == '1H')
Chr2<-subset(genmap_northam, genmap_northam $chromosome == '2H')
Chr3<-subset(genmap_northam, genmap_northam $chromosome == '3H')
Chr4<-subset(genmap_northam, genmap_northam $chromosome == '4H')
Chr5<-subset(genmap_northam, genmap_northam $chromosome == '5H')
Chr6<-subset(genmap_northam, genmap_northam $chromosome == '6H')
Chr7<-subset(genmap_northam, genmap_northam $chromosome == '7H')

DOUBLED_DATA_or_chr1<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr1$SNP)]
DOUBLED_DATA_or_chr2<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr2$SNP)]
DOUBLED_DATA_or_chr3<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr3$SNP)]
DOUBLED_DATA_or_chr4<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr4$SNP)]
DOUBLED_DATA_or_chr5<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr5$SNP)]
DOUBLED_DATA_or_chr6<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr6$SNP)]
DOUBLED_DATA_or_chr7<-DOUBLED_DATA_or[,(colnames(DOUBLED_DATA_or) %in% Chr7$SNP)]

dir.create("~/prephase_files")
for (i in 1:7) {
	INPUT<-get(paste("DOUBLED_DATA_or_chr",i,sep=""))
	OUTPUT<-paste("~/prephase_files/NorthAm_prephase_chr_",i,".txt",sep="")

write.table(INPUT,OUTPUT,quote=F,row.names=T,col.names=T,sep="\t")
}

##Now take the output and run Step2_New_make_fasta_NorthAme.pl to create input for fastPhase
# $ Step2_New_make_fasta_NorthAme.pl NorthAm_prephase_chr_1.txt >./INPUT_NorthAm_chr1.txt

#Add header to INPUT_NorthAm_chr
#Add Sample number, SNP number, a row of P SSSSSS , and a row of SNP positions
#number of SNPs in chromosome
n<-(dim(Chr1)[2])
ADDITIONAL<-rbind(c(1:n),rep("S",n))

write.table(ADDITIONAL,"~/Desktop/Temp_head.txt",quote=F,row.names=F,col.names=F,sep="\t")

#Now run fastPhase using the complete INPUT_NorthAm_chr . One chromsome at the time.