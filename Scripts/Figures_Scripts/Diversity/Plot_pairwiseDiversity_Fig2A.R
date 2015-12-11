#Title           :Plot_pairwiseDiversity_Fig2A.R 
#Description     :Plot genome-wide diversity as in Figure 2A
#Author		 	 :A. Poets 
#========================================================================================

rm(list=ls())

##INPUT FILES: BCAP_FST_Stats.txt is a summary table diversity at each population. The values come from the FIS pipeline and Per_SNP_PairwiseDiversity.R

#the table has the following structure
#	BP RT Nsam SegSites AvgPairwiseDiversity       FIS
#	AB  2  239     2103             0.266012 0.9826380
#	AB  6  142     1791             0.257060 0.9835682
#	BA  2  172     2029             0.203045 0.9969652

TABLE<-read.table("BCAP_FST_Stats.txt",header=T)

Table_ready<-cbind(paste(TABLE[,1],TABLE[,2],sep=""),as.data.frame(TABLE[,c(3:6)]))

colnames(Table_ready)<-c("BP","nsam","segSites","AvgPairwiseDiv","FIS")
row.names(Table_ready)<-c("Idaho (2)","Idaho (6)","Busch Ag. (2)","Busch Ag. (6)","Busch Ag. (Int.)","Minnesota (6)","Montana (2)","North Dakota (2)","North Dakota (6)","Oregon (2)","Oregon (6)","Utah (2)","Utah (6)","Virgina (6)","Washington (2)","Washington (6)")

#Order by diversity
Table_ready_or<-Table_ready[order(Table_ready$AvgPairwiseDiv),]

pdf("~/Dropbox/Mohsen_Fst/Figures/Average_pairwiseDiv.pdf",width=9,height=7)
par(mar=c(7,9,2,2))
barplot(t(as.matrix(Table_ready_or[,4])),horiz=T,names.arg=row.names(Table_ready_or),las=1,xlab="Average Pairwise Diversity")
dev.off()