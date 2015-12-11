#Title           :Plot_Manhattan_Fst_GH_RT_BP_Fig3.R 
#Description     :Make manhattan plots for Fst values as in Figure 3
#Author		 	 :A. Poets 
#========================================================================================

rm(list=ls())

##INPUT FILES: Genetic map and annotations are in Datasets directory.
#Add genetic map information
genmap<-read.table("GeneticMap_T3_020315",header=T)
#Annotations for all BOPA SNPs
ANNOTATIONS<-read.csv("Barley_Annotations.txt",sep="\t")

#Files are the resutls from Fst_six_two_row.R, Fst_winter_spring.r, Fst_BreedingPopulation.R
FST1<-read.table("Fst_winter_spring.txt",header=F)
FST2<-read.table("Fst_RT.txt",header=F)
FST3<-read.table("Fst_BP.txt",header=F)

MAIN<-c("Growth Habit","Row Type", "Breeding Program")
FILE<-c("Fst_GH","Fst_RT","Fst_BP")



for (i in 1:3){

hfst<-get(paste("FST",i,sep=""))
dim(hfst)
head(hfst)

#order SNPs by NAME
hfst_or<-hfst[order(hfst[,1]),]
colnames(hfst_or)<-c("SNPs","FST")

#sort genmap by SNP name
genmap_or<-genmap[order(genmap$SNP),]
#get SNPs that were used in the Fst calculation
genmap_fst<-genmap_or[(genmap_or$SNP %in% hfst_or[,1]),]

if (identical(as.character(genmap_fst$SNP), as.character(hfst_or$SNPs)) == FALSE ) stop (print ("ERROR different SNPs involved"))

#Combine map and FST results AND order by genetic position using the Index column
Fst_order_snp<-cbind(as.data.frame(genmap_fst),as.data.frame(hfst_or))
Fst_order_snp_or<-Fst_order_snp[order(Fst_order_snp$Index),]


chr1<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome =='1H')
chr2<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome=='2H')
chr3<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome=='3H')
chr4<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome=='4H')
chr5<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome=='5H')
chr6<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome=='6H')
chr7<-subset(Fst_order_snp_or, Fst_order_snp_or$chromosome=='7H')

## PLOT

pdf(paste("Fst_manhattanplots/",FILE[i],".pdf",sep=""), width=8, height=5)
par(mar=c(4.5,4,9,4))
plot(Fst_order_snp_or$Cumulative, Fst_order_snp_or$FST,cex=0.45,xlab="Linkage Group",xaxt='n', ylab="Fst", main= paste("Fst ",MAIN[i],sep=""))
points(chr1$Cumulative,chr1$FST, cex=0.45,col='red')
points(chr2$Cumulative,chr2$FST,cex=0.45,col='purple')
points(chr3$Cumulative,chr3$FST,cex=0.45,col='green')
points(chr4$Cumulative,chr4$FST,cex=0.45,col='blue')
points(chr5$Cumulative,chr5$FST,cex=0.45,col='yellow')
points(chr6$Cumulative,chr6$FST,cex=0.45,col='dark green')
points(chr7$Cumulative,chr7$FST,cex=0.45,col='orange')

#To plot the SNPs that fall in a known gene
axis(side=1,at=c(median(chr1$Cumulative),median(chr2$Cumulative),median(chr3$Cumulative),median(chr4$Cumulative),median(chr5$Cumulative),median(chr6$Cumulative),median(chr7$Cumulative)), labels=c("1H","2H","3H","4H","5H","6H","7H"),font=3)

#Draw a line at 97.5%quantile
abline(h=(quantile(Fst_order_snp_or$FST,probs=0.975,na.rm=T)),col="red",lty=4)

dev.off()


### 2. Make table with annotations for SNPs above 97.5% quantile
threshold<-quantile(Fst_order_snp_or$FST,probs=0.975,na.rm=T)
Snp_sig<-subset(hfst_or, hfst_or[,2] >threshold)

SNP_annotation<-ANNOTATIONS[(ANNOTATIONS$SNPName %in% Snp_sig$SNPs),]

#Grab important columns SNPname , Position, and Silent to print tables
SNP_annotation_table<-SNP_annotation[,c("SNPName","GeneShortName","Position","Silent")]

#Order SNPs based on genetic map
SNP_anot_map<-genmap_fst[(genmap_fst$SNP %in% SNP_annotation_table$SNPName),]
if (identical(as.character(SNP_anot_map[,1]),as.character(SNP_annotation_table[,1])) == FALSE) stop (print("SNPs are in different order.Cannot be combined with annotations"))

SNP_annotation_table_genmap<-cbind(as.data.frame(SNP_anot_map[,c(1,3,2,5)]),as.data.frame(SNP_annotation_table[,c(2,3,4)]))

#Add Fst values
if (identical(as.character(Snp_sig[,1]),as.character(SNP_annotation_table_genmap[,1])) == FALSE)stop(print("SNPs are different between significant and annotations"))

SNP_annotation_table_genmap_fst<-cbind(SNP_annotation_table_genmap,as.data.frame(Snp_sig[,2]))
#Sort by chromosome/position
SNP_annotation_table_genmap_or<-SNP_annotation_table_genmap_fst[order(SNP_annotation_table_genmap_fst $Index),]
#Remove INDEX
SNP_annotation_table_genmap_or<-SNP_annotation_table_genmap_or[,-c(4)]

colnames(SNP_annotation_table_genmap_or)<-c("SNPname","Linkage_Group","Position_(cM)","Gene_Short_Name","Position_in_gene","Silent","FST")

SNP_annotation_table_genmap_or[,1]<-sub('X12_',"12_", SNP_annotation_table_genmap_or[,1])
SNP_annotation_table_genmap_or[,1]<-sub('X11_',"11_", SNP_annotation_table_genmap_or[,1])

write.table(SNP_annotation_table_genmap_or,paste(FILE[i],"_sig.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}