rm(list=ls())

##INPUT FILES: from Datasets directory and Phased genotypes
genmap<-read.table("GeneticMap_T3_020315",header=T)
head(genmap)

for (i in 1:7){
	
	#Phased genotypes after the Phasing pipeline
	NORTHAM <-read.table(paste("NorthAm_phased_chr",i,".txt",sep=""),header=T,row.names=1)
	
	##FILE FORMATTING:
	replaceNA<-function(dat){
		dat[is.na(dat)]<-'0'
		return(dat)
		}
	genotypes<-as.data.frame(apply(NORTHAM,2,replaceNA))
	
	#Separate the two haplotypes
	
	hap1<- seq(1,nrow(genotypes),by=2)
	haplotype1<- genotypes[hap1,]
	
	hap2<- seq(2,nrow(genotypes),by=2)
	haplotype2<- genotypes[hap2,]
	
	#combine haplotypes to appear SNP1_allele1 SNP1_allele2
	
	removed<-cbind(haplotype1,haplotype2)
	removed[1:8,1:10]
	dim(removed)
	removed_or<-removed[,order(colnames(removed))]
	dim(removed_or)
	removed_or[1:8,1:10]
	
	#Sort SNPs by genetic position:Turn the table, sort the SNPs based on genetic map, 
	#and turn the table again to have SNPs in columns and Samples in rows
	t_removed_or<-as.data.frame(t(removed_or))
	
	
	#select SNPs present in genetic map and in genotyping data
	
	genmap_part<-genmap[(genmap$SNP %in% row.names(t_removed_or)),]
	
	Index_hap_odd<-seq(1,((dim(NORTHAM)[2])*2),2)
	genmap_odd<-cbind(genmap_part,Index_hap_odd)
	colnames(genmap_odd)<-c("SNP"  , "cm", "Chro", "Cumulative", "index","Index_hap")
	
	Index_hap_even<-seq(2,((dim(NORTHAM)[2])*2),2)
	genmap_even<-cbind(paste(genmap_part[,1],".1",sep=""),as.data.frame(genmap_part[,2:5]),Index_hap_even)
	colnames(genmap_even)<-c("SNP"  , "cm", "Chro", "Cumulative", "index","Index_hap")
	
	genmap_all<-rbind(genmap_odd,genmap_even)
	genmap_all_in<-genmap_all[order(genmap_all$Index_hap),]
	
	rownames(genmap_all_in)<-genmap_all_in$SNP
	
	genmap_all_in_sort<-genmap_all_in[order(row.names(genmap_all_in)),]
	if (identical(rownames(t_removed_or),rownames(genmap_all_in_sort)) == FALSE) stop (print("Error! SNPs are in different order"))
	genotypes_SNP<-cbind(as.data.frame(genmap_all_in_sort), as.data.frame(t_removed_or))
	rownames(genotypes_SNP)<-row.names(t_removed_or)
	
	#sort SNPs by Index_haplotype, this column has the right genetic order for both haplotypes
	genotypes_SNPor<-genotypes_SNP[order(genotypes_SNP$Index_hap),]
	
	#Remove extra columns for SNP ordering and turn the table to have SNPs in columns
	genotypes_SNPor<-genotypes_SNPor[,-c(1:6)]
	genotypes_SNPor_t<-as.data.frame(t(genotypes_SNPor))
	
	#Create extra columns for PED file
	ZERO<-rep("0",(dim(genotypes_SNPor_t)[1]))
	NINE<-rep("-9",(dim(genotypes_SNPor_t)[1]))
	
	PED<-cbind(as.data.frame(ZERO),as.data.frame(row.names(genotypes_SNPor_t)),as.data.frame(ZERO),as.data.frame(ZERO),as.data.frame(ZERO),as.data.frame(NINE),as.data.frame(genotypes_SNPor_t))
	
	#write out input file for Plink for each linkage group
	OUTPUT<-paste("NorthAm_chr",i,".ped",sep="")
	write.table(PED,OUTPUT,quote=F,row.names=F,col.names=F,sep="\t")
	
	
	#Create a.map file for Plink for each linkage group
	snp_single<-cbind(as.data.frame(genmap_part[,c(3,1,2)]),round(genmap_part[,2],digits=0))
	head(snp_single)
	
	#Replace "1H" for just "1" (remove H)
	snp_single[,1]<-gsub(paste(i,"H",sep=""),i, snp_single[,1])
	
	OUTPUT<-paste("chr",i,".map",sep="")
	write.table(snp_single,OUTPUT,quote=F,row.names=F,col.names=F,sep="\t")
}


#Now use Plink_IBS_input_step2.R to separate the files in SNP windows for plink