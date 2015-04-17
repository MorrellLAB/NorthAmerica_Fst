#Read results from Hierarchical Fst. Comparisons winter vs spring. Plot Manhattan plot
rm(list=ls())
#For GH and RT and BP FST
hfst<-read.table("~/Documents/Anita/mohsen_fst/Analysis/for_Fst/Focal_spring_winter/All_winter_all_spring/Fst_winter_spring_ALL.txt",header=F)

dim(hfst)
head(hfst)

#order SNPs by NAME
hfst_or<-hfst[order(hfst[,1]),]
colnames(hfst_or)<-c("SNPs","FST")
#Add genetic map information
genmap<-read.table("~/Documents/Anita/mohsen_fst/Data/Original_QC_diploid_AB/GeneticMap_T3_020315",header=T)
head(genmap)
#sort genmap by SNP name
genmap_or<-genmap[order(genmap$SNP),]
#get SNPs that were used in the Fst calculation
genmap_fst<-genmap_or[(genmap_or$SNP %in% hfst_or[,1]),]

if (identical(dim(genmap_fst)[1], dim(hfst_or)[1])) {print ("Good! Genetic map and SNPs in Fst match")} else {print ("ERROR different SNPs involved")}

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

pdf("~/Documents/Anita/Mohsen_Fst/Plots/Fst_GH.pdf", width=8, height=5)
par(mar=c(4.5,4,9,4))
plot(Fst_order_snp_or$Cumulative, Fst_order_snp_or$FST,cex=0.45,xlab="",xaxt='n', ylab="Fst", main="Fst Growth Habit")
points(chr1$Cumulative,chr1$FST, cex=0.45,col='red')
points(chr2$Cumulative,chr2$FST,cex=0.45,col='purple')
points(chr3$Cumulative,chr3$FST,cex=0.45,col='green')
points(chr4$Cumulative,chr4$FST,cex=0.45,col='blue')
points(chr5$Cumulative,chr5$FST,cex=0.45,col='yellow')
points(chr6$Cumulative,chr6$FST,cex=0.45,col='dark green')
points(chr7$Cumulative,chr7$FST,cex=0.45,col='orange')

#To plot the SNPs that fall in a known gene


axis(side=1,at=c(55.30,218.3,392.1,545.5, 703.2,863.6, 1022.0), labels=c("1H","2H","3H","4H","5H","6H","7H"),font=3)
#Draw a line at 97.5%quantile
abline(h=(quantile(Fst_order_snp_or$FST,probs=0.975,na.rm=T)),col="red",lty=4)

if(FALSE){
#Centromere positions according to Muñoz-Amatriaín et al 2011
#abline(v=c(212.64,391.97,546.17,661.4,864.37,1024.84),col="gray")

# Important genes from file New_coordinates_for Manhattan pots Fst GH and RT
#axis for growth habit related genes
axis(3,at=c(21.3,92.6,184.42,616.93,705.53,734.88,982.12),labels=c("Fr-H3","Ppd-H2","Ppd-H1","Vrn-H2","Fr-H2","Fr-H1/ Vrn-H1","Vrn-H3"), las=3,cex.axis=0.55, font=3)

#axis for Row Type related genes
#axis(3,at=c(236.81,515.12),labels=c("Vrs1","Vrs5"), las=3,cex.axis=0.55,font=3)


#Axis for focal Fst
axis(3,at=c(21.3,92.6,184.42,236.81,515.12,616.93,705.53,734.88,982.12),labels=c("Fr-H3","Ppd-H2","Ppd-H1","Vrs1","Vrs5","Vrn-H2","Fr-H2","Fr-H1/ Vrn-H1","Vrn-H3"), las=3,cex.axis=0.55,font=3)

#Select SNPs in genes for growth habit
SNPsingenes<-c("X12_30930","X12_30239","X12_30239","X12_30873","X12_30894","X12_30871","X12_31163","X12_20187","X11_20653","X11_11490","X12_20403","X11_21314","X11_20749")
SNPsIngenes_data<-Fst_order_snp_or[(Fst_order_snp_or$SNP %in% SNPsingenes),]
#which of these genes are above the 0.975% tail?
cutoff<-quantile(Fst_order_snp_or$FST,probs=0.975,na.rm=T)

passed_cutoff<-subset(SNPsIngenes_data, SNPsIngenes_data$FST >= cutoff)
points(passed_cutoff $Cumulative, passed_cutoff$FST, pch=19,cex=0.45,col='black')
}

dev.off()
