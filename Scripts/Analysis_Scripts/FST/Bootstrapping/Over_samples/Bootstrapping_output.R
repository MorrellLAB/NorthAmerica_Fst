

#Title: Bootstrapping_output.R
#Description: Get mean Fst values after 100 iterations.
#             Plot Fst distribution for one iteration and for 100 iterations.
#Author: Ana Poets
#Note: Output file bootstrapping by Breeding Populations, growth and row types.

rm(list=ls())

##INPUT FILES:
#Add genetic map information
genmap<-read.table("GeneticMap_T3_020315",header=T)
		
#Annotations for all BOPA SNPs
ANNOTATIONS<-read.csv("Barley_Annotations.txt",sep="\t")
	
#Files Fst values from bootstrapping iterations
COMPARISON_boot<-c("Fst_BP_","Fst_RT_","Fst_winter_spring_")

#Files Fst outliers from single iteration.
OUT_sing<-c("Fst_BP_sig.txt","Fst_RT_sig.txt","Fst_GH_All_sig.txt")

#Fst results using entire data set, as in Figure 3.
	One_BP<-read.table("~/All_data/Fst_BP.txt")
	One_RT<-read.table("~/All_data/Fst_RT.txt")
	One_GH<-read.table("~/All_data/Fst_winter_spring_ALL.txt")
	

#Header for plots

MAIN<-c("Breeding populations","Spike type", "Growth habit")

#Results into variables (names)
VAR<-c("BP","RT","GH")


##PROCESSING DATA:
#Table to count number of significant SNPs in a single iteration, and in average across 100 iterations
Table_summary_fst<-matrix(NA,nrow=3,ncol=4)
colnames(Table_summary_fst)<-c("Outliers_single","Outliers_100","Shared","mean variance")
row.names(Table_summary_fst)<-c("Breeding populations","Spike type", "Growth habit")

#Number of iterations
N_ITER<-100

for (c in 1:length(COMPARISON_boot)){
	DATA <-read.table(paste(COMPARISON_boot[c],1,".txt",sep=""))
	FST_results<-matrix(NA,ncol=N_ITER,nrow=(dim(DATA)[1]))
	row.names(FST_results)<-DATA[,1]
	for (b in 1:N_ITER){
		DATA2 <-read.table(paste(COMPARISON_boot[c],b,".txt",sep=""))
		FST_results[,b]<-DATA2[,2]
	}
	
	#Calculate mean and variance in 100 iterations
	average_fst<-as.data.frame(apply(FST_results,1,mean))
	VARIANCE<-function(dat){
	variance<-var(as.numeric(dat),na.rm=T)
	return(variance)
	}

	fst_VAR<-apply(FST_results,1, VARIANCE)
	fst_summary_ite<-cbind(average_fst, fst_VAR)
	
	
	write.table(average_fst,paste("~/averages/Fst_ave_",MAIN[c],".txt",sep=""),quote=F,row.names=T,col.names=F,sep="\t")
	
	#merge the fst values for all the iterations 
	m.data<-as.vector(FST_results)
   assign(paste("m_data_",VAR[c],sep=""),m.data)
	
##PLOTTING:
#make a manhattan plot and a table with significant outliers using mean values
	
		
				
		hfst<-fst_summary_ite
		dim(hfst)
		head(hfst)
		
		#order SNPs by NAME
		hfst_or<-hfst[order(row.names(hfst)),]
		colnames(hfst_or)<-c("FSTmean","FSTvar")
		
		#sort genmap by SNP name
		genmap_or<-genmap[order(genmap$SNP),]
		#get SNPs that were used in the Fst calculation
		genmap_fst<-genmap_or[(genmap_or$SNP %in% row.names(hfst_or)),]
		
		if (identical(as.character(genmap_fst$SNP), as.character(row.names(hfst_or))) == FALSE ) stop (print ("ERROR different SNPs involved"))
		
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
		
		##Plot mean Fst values. The size of the points varies by 100X the variation across 100 iterations
		pdf(paste(VAR[c],"mean_var_resampling",".pdf",sep=""), width=8, height=5)
		par(mar=c(4.5,4,9,4))
		plot(Fst_order_snp_or$Cumulative, Fst_order_snp_or$FSTmean,cex=0.45,xlab="Linkage Group",xaxt='n', ylab="Fst",col="white")
		
		for (s in 1:(dim(chr1)[1])){
		points(chr1$Cumulative[s],chr1$FSTmean[s], cex=((chr1$FSTvar[s])*100),col='red')}
		for (s in 1:(dim(chr2)[1])){
		points(chr2$Cumulative[s],chr2$FSTmean[s], cex=((chr2$FSTvar[s])*100),col='purple')}
		for (s in 1:(dim(chr3)[1])){
		points(chr3$Cumulative[s],chr3$FSTmean[s], cex=((chr3$FSTvar[s])*100),col='green')}
		for (s in 1:(dim(chr4)[1])){
		points(chr4$Cumulative[s],chr4$FSTmean[s], cex=((chr4$FSTvar[s])*100),col='blue')}
		for (s in 1:(dim(chr5)[1])){
		points(chr5$Cumulative[s],chr5$FSTmean[s], cex=((chr5$FSTvar[s])*100),col='yellow')}
		for (s in 1:(dim(chr6)[1])){
		points(chr6$Cumulative[s],chr6$FSTmean[s], cex=((chr6$FSTvar[s])*100),col='dark green')}
		for (s in 1:(dim(chr7)[1])){
		points(chr7$Cumulative[s],chr7$FSTmean[s], cex=((chr7$FSTvar[s])*100),col='orange')}
		
		#To plot the SNPs that fall in a known gene
		axis(side=1,at=c(median(chr1$Cumulative),median(chr2$Cumulative),median(chr3$Cumulative),median(chr4$Cumulative),median(chr5$Cumulative),median(chr6$Cumulative),median(chr7$Cumulative)), labels=c("1H","2H","3H","4H","5H","6H","7H"),font=3)
		
		#Draw a line at 97.5%quantile
		abline(h=(quantile(Fst_order_snp_or$FSTmean,probs=0.975,na.rm=T)),col="red",lty=4)
		
		#Centromere position according to Table 6 in Munoz et al 2011.
		#abline(v=c(210,382.76,552.49,685.17,893,1055))
		dev.off()
		
		
		### 2. Make table with annotations for SNPs above 97.5% quantile
		threshold<-quantile(Fst_order_snp_or$FSTmean,probs=0.975,na.rm=T)
		Snp_sig<-subset(hfst_or, hfst_or[,1] >threshold)
		
		SNP_annotation<-ANNOTATIONS[(ANNOTATIONS$SNPName %in% row.names(Snp_sig)),]
		
		#Grab important columns SNPname , Position, and Silent to print tables
		SNP_annotation_table<-SNP_annotation[,c("SNPName","GeneShortName","Position","Silent")]
		
		#Order SNPs based on genetic map
		SNP_anot_map<-genmap_fst[(genmap_fst$SNP %in% SNP_annotation_table$SNPName),]
		if (identical(as.character(SNP_anot_map[,1]),as.character(SNP_annotation_table[,1])) == FALSE) stop (print("SNPs are in different order.Cannot be combined with annotations"))
		
		SNP_annotation_table_genmap<-cbind(as.data.frame(SNP_anot_map[,c(1,3,2,5)]),as.data.frame(SNP_annotation_table[,c(2,3,4)]))
		
		#Add Fst values
		if (identical(row.names(Snp_sig),as.character(SNP_annotation_table_genmap[,1])) == FALSE)stop(print("SNPs are different between significant and annotations"))
		
		SNP_annotation_table_genmap_fst<-cbind(SNP_annotation_table_genmap,as.data.frame(Snp_sig[,1]))
		#Sort by chromosome/position
		SNP_annotation_table_genmap_or<-SNP_annotation_table_genmap_fst[order(SNP_annotation_table_genmap_fst $Index),]
		#Remove INDEX
		SNP_annotation_table_genmap_or<-SNP_annotation_table_genmap_or[,-c(4)]
		
		colnames(SNP_annotation_table_genmap_or)<-c("SNPname","Linkage_Group","Position_(cM)","Gene_Short_Name","Position_in_gene","Silent","FST")
		
		SNP_annotation_table_genmap_or[,1]<-sub('X12_',"12_", SNP_annotation_table_genmap_or[,1])
		SNP_annotation_table_genmap_or[,1]<-sub('X11_',"11_", SNP_annotation_table_genmap_or[,1])
		
		write.table(SNP_annotation_table_genmap_or,paste("Fst_",VAR[c],"_mean_resampling", "_sig.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
		
		SIG_SNPs_single<-read.table(paste(OUT_sing[c],sep=""),header=T)
		
		#Count how many of SNPs identified as outliers in one iteration are outliers after 100 iterations.
		Table_summary_fst[c,1]<-(dim(SIG_SNPs_single)[1])
		Table_summary_fst[c,2]<-dim(SNP_annotation_table_genmap_or)[1]
		Table_summary_fst[c,3] <-length(intersect(SNP_annotation_table_genmap_or[,1],SIG_SNPs_single[,1]))
		Table_summary_fst[c,4]<-mean(Snp_sig$FSTvar)
}	
write.table(Table_summary_fst,"Table_summary_fst.txt",quote=F,row.names=T,col.names=T,sep="\t")

	