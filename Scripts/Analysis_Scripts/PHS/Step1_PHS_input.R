#Input for Shichen's PHS, analysis based on genome-wide haplotype sharing levels.
rm(list=ls())

#INPUT FILES: In Datasets directory unless indcated otherwise

#Genotype data after running the Phasing pipeline. Use phased heterozygotes, NA for missing data.
#Put all chromosomes together as we are going to do a Genome-wide PHS analysis
CHR1<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr1.txt",header=T,row.names=1)
CHR2<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr2.txt",header=T,row.names=1)
CHR3<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr3.txt",header=T,row.names=1)
CHR4<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr4.txt",header=T,row.names=1)
CHR5<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr5.txt",header=T,row.names=1)
CHR6<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr6.txt",header=T,row.names=1)
CHR7<-read.table("~/final_Phased_diploid_wNA/NorthAm_phased_chr7.txt",header=T,row.names=1)

NorthAm<-cbind(as.data.frame(CHR1),as.data.frame(CHR2),as.data.frame(CHR3),as.data.frame(CHR4),as.data.frame(CHR5),as.data.frame(CHR6),as.data.frame(CHR7))


##Divide by 16 populations

#First we need to double the Samples information file
Samples_info<-read.table("samples_information.txt",header=T)

Samples_info_2<-cbind(paste(Samples_info$Alias,"_2",sep=""),as.data.frame(Samples_info[,2:5]))
colnames(Samples_info_2)<-colnames(Samples_info)

#Put samples info togheter 
Samples_info_doubble<-rbind(Samples_info,Samples_info_2)

#Get sample info only for samples that passed the QC

Sample_info_pass<-Samples_info_doubble[(Samples_info_doubble$Alias %in% row.names(NorthAm)),]


#By breeding program
AB2<-subset(Sample_info_pass, Sample_info_pass $Program == 'AB'& Sample_info_pass $Row_type == '2')
AB6<-subset(Sample_info_pass, Sample_info_pass $Program == 'AB'& Sample_info_pass $Row_type == '6')

BA2<-subset(Sample_info_pass, Sample_info_pass $Program == 'BA'& Sample_info_pass $Row_type == '2')
BA6<-subset(Sample_info_pass, Sample_info_pass $Program == 'BA'& Sample_info_pass $Row_type == '6')

BAI2<-subset(Sample_info_pass, Sample_info_pass $Program == 'BAI'& Sample_info_pass $Row_type == '2')

MT2<-subset(Sample_info_pass, Sample_info_pass $Program == 'MT'& Sample_info_pass $Row_type == '2')

MN6<-subset(Sample_info_pass, Sample_info_pass $Program == 'MN'& Sample_info_pass $Row_type == '6')

N2<-subset(Sample_info_pass, Sample_info_pass $Program == 'N2'& Sample_info_pass $Row_type == '2')
N6<-subset(Sample_info_pass, Sample_info_pass $Program == 'N6'& Sample_info_pass $Row_type == '6')

OR2<-subset(Sample_info_pass, Sample_info_pass $Program == 'OR'& Sample_info_pass $Row_type == '2')
OR6<-subset(Sample_info_pass, Sample_info_pass $Program == 'OR'& Sample_info_pass $Row_type == '6')

UT2<-subset(Sample_info_pass, Sample_info_pass $Program == 'UT'& Sample_info_pass $Row_type == '2')
UT6<-subset(Sample_info_pass, Sample_info_pass $Program == 'UT'& Sample_info_pass $Row_type == '6')

VT6<-subset(Sample_info_pass, Sample_info_pass $Program == 'VT'& Sample_info_pass $Row_type == '6')

WA2<-subset(Sample_info_pass, Sample_info_pass $Program == 'WA'& Sample_info_pass $Row_type == '2')
WA6<-subset(Sample_info_pass, Sample_info_pass $Program == 'WA'& Sample_info_pass $Row_type == '6')


#Get genotyping data for each partition

PARTITIONS<-c("AB2", "AB6", "BA2", "BA6", "BAI2","MT2","MN6","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

for (p in 1:length(PARTITIONS)){
	PART<-get(PARTITIONS[p])
	
	GENOTYPE<-NorthAm[(row.names(NorthAm) %in% PART$Alias),]
	
	write.table(GENOTYPE,paste(PARTITIONS[p],"_genotypes.txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
	
}

#Get Genetic map for all the SNPs in the phased dataset

genmap<-read.table("GeneticMap_T3_020315",header=T)

NorthAm_genmap<-genmap[(genmap$SNP %in% colnames(NorthAm)),]

#Order columns: chr_id SNP_name and genetic_distance (no cummulative)

NorthAm_genmap_or<-as.data.frame(NorthAm_genmap[,c(3,1,2)])

colnames(NorthAm_genmap_or)<-c("chr_id","SNP_name","genetic_distance")

write.table(NorthAm_genmap_or,"All_genetic_distances.txt",quote=F,row.names=F,col.names=T,sep="\t")

#Collect Cumulative cM position to be used while plotting
NorthAm_genmap_or2<-as.data.frame(NorthAm_genmap[,c(3,1,2,4)])
colnames(NorthAm_genmap_or2)<-c("chr_id","SNP_name","genetic_distance","Accumulative_cm")
write.table(NorthAm_genmap_or,"All_genetic_distances_forPlot.txt",quote=F,row.names=F,col.names=T,sep="\t")

#Run calculate_PHS_GW.pl to calculate PHS for each partition


