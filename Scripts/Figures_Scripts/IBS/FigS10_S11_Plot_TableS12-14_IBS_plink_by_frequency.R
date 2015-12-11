#Title           :FigS10_S11_Plot_TableS12-14_IBS_plink_by_frequency.R 
#Description     :Plot IBS from Plink output as in Figures S10 and S11.
#				  Make tables of IBS frequency (Tables S12-S14)
#Author		 	 :A. Poets 
#Note		 	 :Input file comes from Step4_Plink_out_table.R	
#========================================================================================

# Use one population as base line. The Y-axis would be the frequency of shared individuals with other Breeding population evaluated. 
#The thick of the bar for the other programs indicate which frequency of individuals of that program share with the base population.

rm(list=ls())

#Number of SNPs used for each Windows size used for IBS
Windows<-100

DATA_freq <-read.table(paste(Windows,"_SNP_Frequency_shared_seg.txt",sep=""),header=T)

#Create a directory to place the figures
dir.create(paste("~/Figures/",Windows,"_SNPs",sep=""))

#Thinkness of bar represent frequency. So, determine how thick each bar should be
	THICKNESS_VALUE<-function(dat){
		if (findInterval(dat, c(0,0.10) )== 1) (dat<-'1') else
		if (findInterval(dat, c(0.10,0.20) )== 1) (dat<-'2') else
		if (findInterval(dat, c(0.20,0.30) )== 1)  (dat<-'3') else
		if (findInterval(dat, c(0.30,0.40) )== 1) (dat<-'4') else
		if (findInterval(dat, c(0.4,0.50) )== 1)  (dat<-'5') else
		if (findInterval(dat, c(0.5,0.6) )== 1)  (dat<-'6') else
		if (findInterval(dat, c(0.6,0.7) )== 1) (dat<-'7') else
		if (findInterval(dat, c(0.7,0.8) )== 1)  (dat<-'8') else
		if (findInterval(dat, c(0.8,0.9) )== 1)  (dat<-'9')else
		if (findInterval(dat, c(0.9,11) )== 1) (dat<-'10') 
		
	}
	
	
#Choose any row that has the base population


BP<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")


#Make a list of colors for each population

COLOR<-data.frame(c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6"),c("blue","cornflowerblue","darkolivegreen3","darkgreen","green","lightgoldenrod4","gray","purple","slateblue","pink","deeppink","cyan","darkturquoise","deepskyblue4","orange","darkorange3"))

colnames(COLOR)<-c("BP","Color")

legend_table<-data.frame(Pops=c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6"),Names=c("Idaho (2)","Idaho (6)","Busch Ag.(2)","Busch Ag.(6)","Busch Ag. (Int.)","Minnesota (6)","Montana (2)","North Dakota (2)","North Dakota (6)","Oregon (2)","Oregon (6)","Utah (2)","Utah (6)","Virginia (6)","Washington (2)","Washington (6)"),Colors=c("blue","cornflowerblue","darkolivegreen3","darkgreen","green","lightgoldenrod4","gray","purple","slateblue","pink","deeppink","cyan","darkturquoise","deepskyblue4","orange","darkorange3"))


SHARED<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")
for (p in 1:length(SHARED)){
	
	POP_BASE<-SHARED[p]
	POP1<-subset(DATA_freq, DATA_freq $Pop1 == POP_BASE)
	POP2_shared<-subset(DATA_freq, DATA_freq $Pop2 == POP_BASE)
	#Order pops in POP2 to have the based population in columns 1 and 3
	POP2_or<-cbind(as.data.frame(POP2_shared[,c(2,1,4,3,5,6,7)]))
	#Change col names to be the same as POP1 so we can join these data sets
	colnames(POP2_or)<-colnames(POP1)
	
	All_pop<-rbind(POP1,POP2_or)
	
	#Select samples to put in the legend
	Populations_used<-unique(sort(All_pop[,4]))
	legend_use<-legend_table[(legend_table[,1] %in% Populations_used),]

	#PLOT
	pdf(paste("~/Figures/",Windows,"_SNPs/",SHARED[p],"_to_all.pdf",sep=""),width=13,height=7)
	par(mar=c(5.5,4.5,4,15),xpd=T)
	
	POPULATION<-which(legend_table[,1] == SHARED[p])
	plot (c(0.31,1108.97),c(0.001,0.001),col="white",cex=0.05,pch=1,ylim=c(0,1),ylab= (paste("Frequency of IBS segments in ",legend_table[POPULATION,2],sep="")),xlab="Linkage Group",xaxt="n")
	
	par(lend=1)
legend(1160,1.04,lwd=c(1:10),c(">0<0.1","0.1<0.2","0.2<0.3","0.3<0.4","0.4<0.5","0.5<0.6","0.6<0.7","0.7<0.8","0.8<0.9","0.9-0.1"),title=expression(bold("Frequency")),cex=0.9,xpd=T)
	
	#Define linkage groups 
	abline(v=c(140,313,496.13,640,830,972),lty=2,xpd=F)
	#Thick marks for each linkage group
	axis(1,at=c(70,226.5,404.075,567.575,735,901,1040.485),label=c("1H","2H","3H","4H","5H","6H","7H"))
	
	#for each pair find which population is interacting, this will help to determine the color to be used
	for (j in 1:dim(All_pop)[1]){
	POP2<-All_pop[j,4]
	COLOR_VAL<-COLOR[(COLOR$BP %in% POP2[1]),]
	
	par(lend=1)
	segments(All_pop[j,6], All_pop[j,1], All_pop[j,7], All_pop[j,1],col=as.character(COLOR_VAL[1,2]),lwd= THICKNESS_VALUE(All_pop[j,2]))
	}
	dev.off()
}
	

#============TABLE==========
#TABLE S12
#Make a table with the number of segments shared with each population/ number of those that are in >50%freq in both pops

BP<-c("AB2","BA2","BAI2","MT2","N2","OR2","UT2","WA2","AB6","BA6","MN6","N6","OR6","UT6","VT6","WA6")

TABLE<-matrix(NA,ncol=length(BP),nrow=length(BP))
colnames(TABLE)<-BP
row.names(TABLE)<-BP

head(DATA_freq)

for (t in 1:length(BP)){
	POP<-subset(DATA_freq, DATA_freq$Pop1 == BP[t])
	totalSegments<-length(unique(sort(POP$Segment)))
	#For each population contributing to the one selected count the number of unique segments contributed, and then count only those segments 
	#contributed more than 40%of donor population to at leat 20% of recipient population
	for (c in 1:length(BP)){
		Contributor<-subset(POP, POP$Pop2 == BP[c])
		
		if (dim(Contributor)>0) {
			sharedSegments<-((length(unique(sort(Contributor$Segment))))/totalSegments)
			sharedSegments <-round(sharedSegments,digits=2)
			
			percentShared40<-subset(Contributor, Contributor$freq_pop1 >=0.2 & Contributor$freq_pop2 >=0.40)
			
			#percentage of  genomic segments shared at the desired proportions
			countSharedHigh<-((length(unique(sort(percentShared40$Segment))))/totalSegments)
			countSharedHigh<-round(countSharedHigh,digits=2)
			resultCount<-paste(sharedSegments, countSharedHigh,sep="/")
			
			#Find position in the table
			ROW<-which(row.names(TABLE) == BP[t])
			COL<-which(colnames(TABLE) == BP[c])
			TABLE[ROW,COL]<-resultCount
		}
	}
	
}


legend_table2<-c("Idaho (2)","Busch Ag.(2)","Busch Ag. (Int.)","Montana (2)","North Dakota (2)","Oregon (2)","Utah (2)","Washington (2)","Idaho (6)","Busch Ag.(6)","Minnesota (6)","North Dakota (6)","Oregon (6)","Utah (6)","Virginia (6)","Washington (6)")
colnames(TABLE)<-legend_table2
row.names(TABLE)<-legend_table2

TABLE<-as.data.frame(TABLE)

write.table(TABLE,paste("~/Documents/Anita/Mohsen_Fst/Tables/TableS_frequencyIBS_",Windows,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")

#=== TABLE S13 & TABLE S14: Make a table of percentage of shared segments "/" average frequency of shared segments

BP<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

TABLE_proportion_freq<-matrix(NA,ncol=length(BP),nrow=length(BP))
colnames(TABLE_proportion_freq)<-BP
row.names(TABLE_proportion_freq)<-BP

head(DATA_freq)

for (t in 1:length(BP)){
	POP<-subset(DATA_freq, DATA_freq$Pop1 == BP[t])
	totalSegments<-length(unique(sort(POP$Segment)))
	#For each population contributing to the one selected count the number of unique segments contributed, and then count only those segments 
	#contributed more than 40%of donor population to at leat 20% of recipient population
	for (c in 1:length(BP)){
		Contributor<-subset(POP, POP$Pop2 == BP[c])
		
		if (dim(Contributor)>0) {
			sharedSegments<-((length(unique(sort(Contributor$Segment))))/totalSegments)
			sharedSegments <-round(sharedSegments,digits=2)
			meanFreq_pop1<-mean(Contributor[,1])
			meanFreq_pop1<-round(meanFreq_pop1,digits=2)
			meanFreq_pop2<-mean(Contributor[,2])
			meanFreq_pop2<-round(meanFreq_pop2,digits=2)
			
			resultCount<-paste(sharedSegments," (",meanFreq_pop1,"/", meanFreq_pop2,")",sep="")
			
			#Find position in the table
			ROW<-which(row.names(TABLE_proportion_freq) == BP[t])
			COL<-which(colnames(TABLE_proportion_freq) == BP[c])
			TABLE_proportion_freq[ROW,COL]<-resultCount
		}
	}
	
}

colnames(TABLE_proportion_freq)<-legend_table[,2]
row.names(TABLE_proportion_freq)<-legend_table[,2]

TABLE_proportion_freq <-as.data.frame(TABLE_proportion_freq)

write.table(TABLE_proportion_freq,paste("TableS_Proportion_frequencyIBS_",Windows,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
