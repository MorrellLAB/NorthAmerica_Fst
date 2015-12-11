#Title           :Step4_Plink_out_table.R
#Description     :Creates a table quantifiying the frequency of shared alleles at a given window between two populations
#Author		 	 :A. Poets
#Date			 :Jun 26, 2015
#Note		     :Requires output files from Step3_Plink_IBS_OUT.sh
#========================================================================================

rm(list=ls())

#Identify the window size used for IBS in PLINK
Windows<-50


##INPUT FILES:unless indicated otherwise the files are in Dataset directory.
#Get the list of segments
	#This tables were created prior to run Plink in Step2.
	segment_list<-read.table(paste(Windows ,"_SNPs_out/List_gral.txt",sep=""))

	#Get segments information
	Seg_info<-read.table(paste("~/Dropbox/Mohsen_Fst/Analysis/for_IBS_PLINK/Input/PLINK/Input/", Windows ,"_SNPs/Segments_info.txt",sep=""))

# Call genetic map
genmap <-read.table("GeneticMap_T3_020315",header=T)

#Select markers used in IBS analysis created in Step2
IBS_map<-read.table(paste(Windows ,"_SNPs/map_cm_", Windows ,"snpALL.map",sep=""),header=F)

genmap_ibs<-genmap[(genmap$SNP %in% IBS_map[,2]),]

#Get file indicating which samples belong to which of the 16 breeding populations
Samples_information<-read.table("samples_information.txt",header=T)
PROGRAMS<-as.data.frame(cbind(paste(Samples_information[,3], Samples_information[,4],sep=""), as.character(Samples_information[,1])))
colnames(PROGRAMS)<-c("PROGRAM","Sample")

#Get the start and end position for each Chrom segment
Seg_info_correct<-matrix(0,ncol=3,nrow=(dim(Seg_info)[1]))
colnames(Seg_info_correct)<-c("Segment","Start","End")
for (c in 1:(dim(Seg_info)[1])){
Seg_info_correct[c,1]<-as.character(Seg_info[c,1])
Seg_info_correct[c,2]<-genmap_ibs[(which(genmap_ibs$SNP == as.character(Seg_info[c,4]))),4]
Seg_info_correct[c,3]<-genmap_ibs[(which(genmap_ibs$SNP == as.character(Seg_info[c,10]))),4]
}
Seg_info_correct <-as.data.frame(Seg_info_correct)

#Make a table with the amount of samples per breeding program
TotalSamples<-as.data.frame(table(PROGRAMS[,1]))

# For each segment test if the individuals are from different populations, if so, assign a color to each individual and make a line with each color

head(segment_list)

##CONSTRUCT TABLE:

BP<-c("AB2","AB6","BA2","BA6","BAI2","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

FREQ_SEQ<-as.data.frame(matrix(ncol=7,nrow=0))
colnames(FREQ_SEQ)<-c("freq_pop1", "freq_pop2", "Pop1", "Pop2", "Segment", "Start", "End")

#for (s in 1:(dim(segment_list)[1])){
IBS_out_all<-read.table(paste(Windows ,"_SNPs/",segment_list[s,1],"_perfectMatch.txt",sep=""))
	
	#Remove combinations with more than 10% missmatch (Col12<0.90)
	More10mistmatch<-which(IBS_out_all[,12] < 0.9)
	if(length(More10mistmatch) > 0) {IBS_out_all_90 <- IBS_out_all[-c(More10mistmatch),]} else {IBS_out_all_90<-IBS_out_all}
	IBS_out<-IBS_out_all_90[,1:4]

	for (l in 1:length(BP)){
	BP_assig<-subset(PROGRAMS, PROGRAMS$PROGRAM == BP[l])
	samples_assigned<-is.element(IBS_out[,2], BP_assig $Sample)
	toChange<-which(samples_assigned ==TRUE)
	if (length(toChange) >0){
	IBS_out[toChange,5]<-c(BP[l])
		}
	samples_assigned<-is.element(IBS_out[,4], BP_assig $Sample)
	toChange4<-which(samples_assigned ==TRUE)
	if (length(toChange4) >0){
	IBS_out[toChange4,6]<-c(BP[l])
		}
	}

	#Remove rows with same breeding program

	SAME_POP<-function(dat){
	if (dat[5] == dat [6]) (CHECK <-c("EQUAL")) else ((CHECK<-c("DIFFERENT")))
	return(CHECK)
	}
	Duplicates<-apply(IBS_out,1,SAME_POP)
	Are_duplicates<-which(Duplicates == "EQUAL")

	#Removing rows duplicates
	IBS_out_unique<-IBS_out[-c(Are_duplicates),]

	#Remove duplicated pairs (pairs of different i.e. BP: A->B == B->A)
	DUPLICATED_OR<-function(dat){
    dat_col<-dat[c(2,4)]
    dat_or<-dat_col[order(dat_col)]
    return(dat_or)
    }
 
	NAMES_ORDERED<-t(as.data.frame(apply(IBS_out_unique, 1, DUPLICATED_OR)))
	
	duplicated_pairs<-which(duplicated(NAMES_ORDERED)== TRUE)

	
	#if exist, remove duplicated pairs
	if (length(duplicated_pairs) >0) {PAIRS_UNIQUE<-IBS_out_unique[-c(duplicated_pairs),]} else IBS_out_unique -> PAIRS_UNIQUE

	#Sometimes, there are (A, B) and (B,A) breediong populations pairs. Since these are basically the same order the BP
	#assignations in alphabitical order for every pair in columns 5 and 6, then select all those pairs to estimate
	#the frequency of individuals in A and in B.

	#Order pairs alphabetically
	DUPLICATED_5_6<-function(dat){
    dat_col<-dat[c(5,6)]
    dat_or<-dat_col[order(dat_col)]
    return(dat_or)
    }
 
	GUIDE_PAIRS<-as.data.frame(t(as.data.frame(apply(PAIRS_UNIQUE, 1, DUPLICATED_5_6))))
	
	#make a list of all unique pairs in the dataset
	list_pairs<-as.data.frame(table(GUIDE_PAIRS))
	#Get pairs that actually exist
	list_pairs_value<-subset(list_pairs, list_pairs$Freq >0)
	
	PAIRS_UNIQUE_or<-cbind(as.data.frame(PAIRS_UNIQUE),as.data.frame(GUIDE_PAIRS))
	colnames(PAIRS_UNIQUE_or)<-c("extra1","sampleBP1","extra2","sampleBP2","BP1","BP2","GUIDE1","GUIDE2")
 	
 	#For each pair shared calculate the frequency of individuals in each population
 	
	 for (sh in 1:(dim(list_pairs_value)[1])) {
	 	
	 	PAIR<-subset(PAIRS_UNIQUE_or ,PAIRS_UNIQUE_or$GUIDE1 == list_pairs_value[sh,1] & PAIRS_UNIQUE_or$GUIDE2 == list_pairs_value[sh,2])
	 	
	 	if (dim(PAIR)[1] >0){
	 	#Put all samples on top of each other and calculate frequencies
	 	pop1<-PAIR[,c(2,5)]
	 	colnames(pop1)<-c("Samples","BP")
	 	pop2<-PAIR[,c(4,6)]
	 	colnames(pop2)<-c("Samples","BP")
			
		Get_frequencies<-rbind(as.data.frame(pop1),as.data.frame(pop2))
		pop1_ind<-subset(Get_frequencies ,Get_frequencies[,2] == list_pairs_value[sh,1])
	 	pop2_ind<-subset(Get_frequencies ,Get_frequencies[,2] == list_pairs_value[sh,2])
	 		
	 	pop1_count<-length(unique(sort(pop1_ind[,1])))
	 	pop2_count<-length(unique(sort(pop2_ind[,1])))
	 	#Get frequencies
	 	freq_pop1<-pop1_count /TotalSamples[which(TotalSamples[,1] == as.character(list_pairs_value[sh,1])),2]
	 	freq_pop2<-pop2_count /TotalSamples[which(TotalSamples[,1] == as.character(list_pairs_value[sh,2])),2]
	 	position_Segment<-which(Seg_info_correct[,1] == segment_list[s,1] )
	 		
	 	freq_segments<-cbind(freq_pop1, freq_pop2,list_pairs_value[sh,1],list_pairs_value[sh,2], Seg_info_correct[position_Segment,])
	 	colnames(freq_segments)<-c("freq_pop1", "freq_pop2", "Pop1", "Pop2", "Segment", "Start", "End")
	 	#Frequency of population to be base in the plot
		FREQ_SEQ<-rbind(FREQ_SEQ, freq_segments)
			}
	 	}
 }
 
 write.table(FREQ_SEQ,paste(Windows,"_SNP_Frequency_shared_seg.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
 