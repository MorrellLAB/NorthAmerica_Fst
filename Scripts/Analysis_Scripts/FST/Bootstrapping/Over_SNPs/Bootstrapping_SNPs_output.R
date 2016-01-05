#Title: Bootstrapping_SNPs_output.R
#Description: Count number of regions that are significant in every iteration that are 
#             also signficant when the whole data set is used to estimate Fst.
#			  Make a manhattan plot for every iteration.
#Author: A.Poets
#==============================================================================

rm(list=ls())

##INPUT FILES:
#Genetic map
genmap<-read.table("GeneticMap_T3_020315",header=T)

#Bootstrapping output files path
BOOTS<-c("All_BP/Fst_BP_","Rowtype_Fst/Fst_RT_","S2_S6_resampling/Fst_GH_")

#Number of regions per partition
N_regions<-c(8,14,16)

#Number of iterations
N_ITERATION<-100

#Directories to put plots
PLOTS<-c("All_BP","Rowtype_Fst","S2_S6_resampling")
MAIN<-c("Breeding populations","Spike type","Groth habit")

PARTITION<-c("BP","RT","GH")
###===Genomic regions start and end position.
 #Regions with significant Fst values in the Breeding populations comparison. Distance in cM. From table S4.
 #Ignoring SNPs in linkage group 1 as there are significant SNPs all across the chromosome.

BP_region1<-c(209.41,217.1)
BP_region2<-c(236.47,239.44)

BP_region3<-c(378.08,384.65)

BP_region4<-c(524.27,525.02)
BP_region5<-c(552.49,556.76)
BP_region6<-c(629.05,641.97)

BP_region7<-c(740.96,742.15)

BP_region8<-c(893.85,895.37)

#Regions with significant Fst values in the Growth habit comparison. Distance in cumulative cM. From table S5.
GH_region1<-c(91.73,93.48)
GH_region2<-c(140.69,142.55)

GH_region3<-c(168.54,175.12)
GH_region4<-c(196.3,227.56)
GH_region5<-c(292.48,316.13)

GH_region6<-c(493.15,494.39)

GH_region7<-c(517.8,518.54)
GH_region8<-c(552.49,553.15)
GH_region9<-c(629.05,634.97)

GH_region10<-c(694.27,702.97)
GH_region11<-c(710.42,727.21)
GH_region12<-c(740.96,756.27)
GH_region13<-c(763.11,788.66)

GH_region14<-c(928.36,928.36)
GH_region15<-c(964.3,964.3)

GH_region16<-c(1048.41,1062.93)

#Regions with significant Fst values in the Row type comparison. Distance in cM. From table S6.
#Ignoring SNPs in linkage group 1 as there are significant SNPs all across the chromosome.

RT_region1<-c(236.47,239.44)

RT_region2<-c(379.39,380.08)
RT_region3<-c(424.84,425.57)
RT_region4<-c(451.57,452.32)

RT_region5<-c(524.27,526.64)
RT_region6<-c(543.24,545.04)
RT_region7<-c(593.32,594.97)
RT_region8<-c(630.82,641.97)

RT_region9<-c(644.67,644.67)
RT_region10<-c(681.54,682.24)
RT_region11<-c(694.27,696.6)
RT_region12<-c(741.61,745.77)

RT_region13<-c(908.69,911.18)

RT_region14<-c(1036.54,1040.85)

#Make a list with all the SNPs that fell between the boundaries for each region in each comparison


for (p in 1:length(PARTITION)){
	#maximum number of regions in any given comparison
	for (r in 1:N_regions[p]){
		Region<-paste(PARTITION[p],"_region",r,sep="")
		if(exists(Region[1]) == FALSE)next
		Boundaries<-get(Region)
		START<-min(which(genmap$Cumulative == Boundaries[1]))
		END<-max(which(genmap$Cumulative == Boundaries[2]))

		REGION_SNPs<-genmap[c(START:END),]
		assign(paste(PARTITION[p],"_regionSNP_",r,sep=""), REGION_SNPs)		
	}	
}

#Take one iteration at the time and count how many REGIONS are significant. One significant SNP is enough to count the region.
#Take one comparison (BP,RT,GH) at the time
for (p in 1:length(PARTITION)){
	TABLE_SUMMARY<-matrix(NA,ncol=3,nrow=N_ITERATION)
	colnames(TABLE_SUMMARY)<-c("Sig. SNPs","Expected No. Regions","Observed No.Regions")
	#take one iteration at the time and identify significant SNPs.
	for (i in 1:N_ITERATION){
		DATA<-read.table(paste(BOOTS[p],i,".txt",sep=""))
		threshold<-quantile(DATA[,2],probs=0.975,na.rm=T)
		SIGNIFICANT_SNPS<-subset(DATA,DATA[,2] >= threshold)
		
		#Look if there are significant SNPs in the expected regions
		COUNT_SIG<-c(0)
		for (R in 1:N_regions[p]){
			expected_region<-paste(PARTITION[p],"_regionSNP_",R,sep="")
			if(exists(expected_region[1]) == FALSE)next
				SNPs_expected<-get(expected_region)
				if (length(intersect(SIGNIFICANT_SNPS[,1], SNPs_expected[,1])) >0) {COUNT_SIG <-COUNT_SIG+1}
			}
		#Number of Recovered Regions	
		TABLE_SUMMARY[i,1]<-dim(SIGNIFICANT_SNPS)[1]
		TABLE_SUMMARY[i,2]<-N_regions[p]
		TABLE_SUMMARY[i,3]<-COUNT_SIG
	}
	AVERAGE<-apply(TABLE_SUMMARY,2,mean)
	Table_final<-rbind(TABLE_SUMMARY,AVERAGE)
	assign(paste(PARTITION[p],"_tableSummary_sig",sep=""), Table_final)

	write.table(get(paste(PARTITION[p],"_tableSummary_sig",sep="")),paste(PARTITION[p],"_TableSumary_RegionsSig.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}


##===Plot every iteration
for (p in 1:length(PARTITION)){
	
	#take one iteration at the time and identify significant SNPs.
	for (i in 1:N_ITERATION){
		hfst <-read.table(paste(BOOTS[p],i,".txt",sep=""))
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
		
		## manhattan plot
		
		pdf(paste(PLOTS[p],"/", PARTITION[p],"-",i,".pdf",sep=""), width=8, height=5)
		par(mar=c(4.5,4,9,4))
		plot(Fst_order_snp_or$Cumulative, Fst_order_snp_or$FST,cex=0.45,xlab="Linkage Group",xaxt='n', ylab="Fst", main= paste("Fst ",MAIN[p],sep=""))
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
	}		
}		

