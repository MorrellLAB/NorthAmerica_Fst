#PLOT PHS RESULTS AFTER RUNNING calculate_PHS_GW.pl 
rm(list=ls())

##INPUT FILE:
Annotations<-read.csv("Barley_Annotations.txt",sep="\t",header=T)
#Table generated from Step1_PHS_input.R
SNP_cumulPosition<-read.table("All_genetic_distances_forPlot.txt",header=T)

##FILE FORMATTING:
#adjusting output from PHS to get ready for plot
breeding_program<-as.data.frame(c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6"))

#For each breeding population

for (i in 1:(dim(breeding_program)[1])){
pop_analyzed<-breeding_program[i,]

#Create a directory where the combine output files per population should go.
dir.create("~/",pop_analyzed,sep="")
##call the results from PHS analysis
chr1<-read.table(paste(pop_analyzed,"/chromosome_1H_phs.txt",sep=""),header=T)
chr2<-read.table(paste(pop_analyzed,"/chromosome_2H_phs.txt",sep=""),header=T)
chr3<-read.table(paste(pop_analyzed,"/chromosome_3H_phs.txt",sep=""),header=T)
chr4<-read.table(paste(pop_analyzed,"/chromosome_4H_phs.txt",sep=""),header=T)
chr5<-read.table(paste(pop_analyzed,"/chromosome_5H_phs.txt",sep=""),header=T)
chr6<-read.table(paste(pop_analyzed,"/chromosome_6H_phs.txt",sep=""),header=T)
chr7<-read.table(paste(pop_analyzed,"/chromosome_7H_phs.txt",sep=""),header=T)

rbind(as.data.frame(chr1),as.data.frame(chr2),as.data.frame(chr3),as.data.frame(chr4),as.data.frame(chr5),as.data.frame(chr6),as.data.frame(chr7)) ->PHS_results

head(PHS_results)

write.table(PHS_results,paste("~/",pop_analyzed, "/",pop_analyzed,"_all_chr_PHS.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
	}


#Call PHS results grouped for all chromosomes in the population
for (i in 1:(dim(breeding_program)[1])){
	pop_analyzed<-breeding_program[i,]
	
	PHS_results<-read.table(paste("~/",pop_analyzed, "/",pop_analyzed,"_all_chr_PHS.txt",sep=""),header=T)
	
	if (identical(as.character(SNP_cumulPosition$SNP_name),as.character(PHS_results$SNP_name)) == FALSE ) stop (print("SNPs are in different order"))
	#Select CHR,SNP_name,cM_cumul, PHS_A,PHS_B,FREQ_A,FREQ_B
	PHS<-cbind(as.data.frame(SNP_cumulPosition [,c(1,2,4)]),as.data.frame(PHS_results[,c(5,6,7,8)]))
	

	#Add Accumulative genetic position to the PHS values. This for plotting purposes.
	Accum_poss_shared<-SNP_cumulPosition[(SNP_cumulPosition $SNP_name %in% PHS_results$SNP_name),]
	if (identical(as.character(PHS_results$SNP_name), as.character(Accum_poss_shared$SNP_name)) == FALSE) stop (print("Error! SNPs are in different order"))
	ACCUMULATIVE<-Accum_poss_shared$cM_cumulative
	cbind(PHS,ACCUMULATIVE)->PHS_IN

	#Select values in the 97.5% quantile for PHS at each allele, and allele freq of >=0.1 (arbitrary)
	cutoff_A<-quantile(PHS $PHS_A,probs=0.975,na.rm=T)
	cutoff_B<-quantile(PHS $PHS_B,probs=0.975,na.rm=T)
	
	significant_A<-subset(PHS, PHS $PHS_A >= cutoff_A & PHS $Freq_A >=0.1)
	significant_B<-subset(PHS, PHS $PHS_B >= cutoff_B & PHS $Freq_B >=0.1)

	#Obtain a list of unique Significan SNPs, when SNP is significant in either Allele_A or Allele_B
	rbind(significant_A,significant_B)->all_sig
	unique(all_sig)->un_all_sig
	colnames(un_all_sig)<-c("Chr","SNP_name","Acum_Genetic_dist","PHS_A","PHS_B","Freq_A","Freq_B")
	unique_significant<-un_all_sig[,c("Chr","SNP_name","Acum_Genetic_dist","PHS_A","PHS_B","Freq_A","Freq_B")]
	
	
	# Get the SNP annotation for significant SNPs
	Annotation_sig<-Annotations[(Annotations$SNPName %in% unique_significant$SNP_name),]
	Annotation_sig_or<-Annotation_sig[order(Annotation_sig$SNPName),]
	unique_significant_or<-unique_significant[order(unique_significant$SNP_name),]
	if (identical(as.character(Annotation_sig_or$SNPName),as.character(unique_significant_or$SNP_name))== FALSE) stop (print ("Error! Annotation and unique SNPs in different order"))
	
	unique_significant_annot<-cbind(as.data.frame(unique_significant_or),as.data.frame(Annotation_sig_or[,c("GeneShortName","Position","Silent","ProductName")]))
	#Sort by genetic position
	unique_significant_annot_or<-unique_significant_annot[order(unique_significant_annot$Acum_Genetic_dist),]
	write.table(unique_significant_annot_or,paste("Significant_",breeding_program[i,],".xls",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
	

	#Calculate which cutoff is lower, to trace the line in the plot
	min_threshold<-min(cutoff_A,cutoff_B)
	
	##GET SIGNIFICANT SNPS AND THEIR BOUNDARIES (Segments)
	PHS_results[(PHS_results$SNP_name %in% unique_significant$SNP_name), ] ->u_sig
	
	head(u_sig)
	
	#Divide the segmenents boundaries for each allele
	u_sig$Block_A1<-lapply(strsplit(as.character(u_sig$Block_A), "\\_") , "[",1)
	u_sig$Block_A2<-lapply(strsplit(as.character(u_sig$Block_A), "\\_") , "[",2)
	
	u_sig$Block_B1<-lapply(strsplit(as.character(u_sig$Block_B), "\\_") , "[",1)
	u_sig$Block_B2<-lapply(strsplit(as.character(u_sig$Block_B), "\\_") , "[",2)
	
	as.data.frame(u_sig) -> u_sig
	
	##Find the difference between the SNPposition in GeneticDistance relative to upper and lower boundaries
	
	#Allele A
	u_sig$Lower_A<-c(u_sig$Genetic_dist) - as.numeric(c(u_sig$Block_A1))
	u_sig$Upper_A<-as.numeric(c(u_sig$Block_A2)) - c(u_sig$Genetic_dist)
	
	#AlleleB
	u_sig$Lower_B<-c(u_sig$Genetic_dist) - as.numeric(c(u_sig$Block_B1))
	u_sig$Upper_B<-as.numeric(c(u_sig$Block_B2)) - c(u_sig$Genetic_dist)
	
	#Select the significant SNPs from PHS output files that have the Accumulative genetic position
	PHS_IN[(PHS_IN$SNP_name %in% u_sig$SNP_name),] ->acc_pos
	ACCUMULATIVE<-acc_pos$ACCUMULATIVE
	if( identical(as.character(acc_pos$SNP_name),as.character(u_sig$SNP_name)) == FALSE) stop (print("Error! SNP names are in different order"))
	cbind(as.data.frame(ACCUMULATIVE),u_sig) ->PHS_accu
	
	#Take the cM position corresponding to the lower and upper boundaries (difference or increase from Accumulative)
	PHS_accu[,1] - PHS_accu$Lower_A -> BlockA1_accum
	PHS_accu[,1] + PHS_accu$Upper_A -> BlockA2_accum
	
	PHS_accu[,1] - PHS_accu$Lower_B -> BlockB1_accum
	PHS_accu[,1] + PHS_accu$Upper_B -> BlockB2_accum
	
	PHS_ready_plot<-cbind(PHS_accu,as.data.frame(BlockA1_accum),as.data.frame(BlockA2_accum),as.data.frame(BlockB1_accum),as.data.frame(BlockB2_accum))

	# Select only the significant points either for allele A or B. Here we want only the cM position in general.
	significant_a<-subset(PHS_ready_plot,PHS_ready_plot$PHS_A >= min_threshold)
	significant_b<-subset(PHS_ready_plot,PHS_ready_plot$PHS_B >= min_threshold)
	
	#Create a population identifier for each allele
	Pop_a<-rep(pop_analyzed,dim(significant_a)[1])
	Pop_b<-rep(pop_analyzed,dim(significant_b)[1])
	
	#Put together the pop_identifier with the significant markers
	
	assign(paste("SIGNIFICANT_",pop_analyzed,"_a",sep=""), cbind(as.data.frame(Pop_a),as.data.frame(significant_a)))
	assign(paste("SIGNIFICANT_",pop_analyzed,"_b",sep=""), cbind(as.data.frame(Pop_b),as.data.frame(significant_b)))
	}
	

	###Put together all the Significant_a for all pops , and another file the significant for b allele.
	
	
	
	All_pops_a <-rbind(as.data.frame(SIGNIFICANT_AB2_a),as.data.frame(SIGNIFICANT_AB6_a),as.data.frame(SIGNIFICANT_BA2_a),as.data.frame(SIGNIFICANT_BA6_a),as.data.frame(SIGNIFICANT_BAI2_a),as.data.frame(SIGNIFICANT_MN6_a),as.data.frame(SIGNIFICANT_MT2_a),as.data.frame(SIGNIFICANT_N2_a),as.data.frame(SIGNIFICANT_N6_a),as.data.frame(SIGNIFICANT_OR2_a),as.data.frame(SIGNIFICANT_OR6_a),as.data.frame(SIGNIFICANT_UT2_a),as.data.frame(SIGNIFICANT_UT6_a),as.data.frame(SIGNIFICANT_VT6_a),as.data.frame(SIGNIFICANT_WA2_a),as.data.frame(SIGNIFICANT_WA6_a))
	
	
	All_pops_b <-rbind(as.data.frame(SIGNIFICANT_AB2_b),as.data.frame(SIGNIFICANT_AB6_b),as.data.frame(SIGNIFICANT_BA2_b),as.data.frame(SIGNIFICANT_BA6_b),as.data.frame(SIGNIFICANT_BAI2_b),as.data.frame(SIGNIFICANT_MN6_b),as.data.frame(SIGNIFICANT_MT2_b),as.data.frame(SIGNIFICANT_N2_b),as.data.frame(SIGNIFICANT_N6_b),as.data.frame(SIGNIFICANT_OR2_b),as.data.frame(SIGNIFICANT_OR6_b),as.data.frame(SIGNIFICANT_UT2_b),as.data.frame(SIGNIFICANT_UT6_b),as.data.frame(SIGNIFICANT_VT6_b),as.data.frame(SIGNIFICANT_WA2_b),as.data.frame(SIGNIFICANT_WA6_b))
	
	#There was some problem with the output looking like a list, so i use the lapply to coarse into a data frame
	All_pops_a_2 <- data.frame(lapply(All_pops_a, as.character), stringsAsFactors=FALSE)
	All_pops_b_2 <- data.frame(lapply(All_pops_b, as.character), stringsAsFactors=FALSE)
	write.table(All_pops_a_2,"Allele_A_all_16_pops.txt",col.names=T,row.names=F,quote=F,sep="\t")
	
	write.table(All_pops_b_2,"Allele_B_all_16_pops.txt",col.names=T,row.names=F,quote=F,sep="\t")
	


#continue in Plot_FigS8_TableS9_PHS_all_pops_Segments.R