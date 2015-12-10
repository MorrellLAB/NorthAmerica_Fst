#After runing Step1_IBS_input_PEDfiles.R 
#THIS CODE WILL:
#1. Remove the SNPs at the end of each chromosome that can not be devided in 20 (20SNP windows)
#2. Create new map files with only the remaining SNPs
#3.	Split each chromosome in XSNP windows for both the new .map and the .ped files.

rm(list=ls())

#1. .map files created in Step1_IBS_input_PEDfiles.R 

#Determine the windows size
SNPcount <-50

#Create a directory for each window size

dir.create(paste(SNPcount,"_SNPs",sep=""))

WORKDIR<-paste(SNPcount,"_SNPs",sep="")

for (i in 1:7){
	##INPUT FILE: output from Step1
	CHR<-read.table(paste("chr",i,".map",sep=""))
	
	CHR_SNP<-floor(dim(CHR)[1]/SNPcount)* SNPcount
	#Write .map files with just the number of SNPs that are devided by 20.
	Chr_map <-as.data.frame(CHR[1:CHR_SNP,])
	Output<-paste(WORKDIR,"/chr",i,"_",SNPcount,"s.map",sep="")
	write.table(Chr_map,Output,,quote=F,row.names=F,col.names=F,sep="\t")
	

	###2. Devide the map and .ped files in segments of 20SNPs
	
	##INPUT FILE:output from Step1
	PED<-read.table(paste("NorthAm_chr",i,".ped",sep=""),header=F)
	
	#Separate sample names from ped file
	
	SAMPLES<-PED[,1:6]
	PED<-PED[,-c(1:6)]
	
	#Select the same number of markers that the new 20s.map has. Remember that the .ped has two columns per marker
	#since we removed the end of each chromosome, we can select only the lenght that reminds times 2, for two alleles per marker
	Trim_ped<-PED[,1:(dim(Chr_map)[1]*2)] 
	dim(Trim_ped)
	
	array_ped<-seq(1, (dim(Trim_ped)[2]),(SNPcount*2))
	array_map<-seq(1, (dim(Chr_map)[1]),(SNPcount))
	
		#Separate map in XSNP map files and .ped files.
		for (s in 1:length(array_map)){
			Map<-Chr_map[(array_map[s]):(array_map[s]+(SNPcount -1)),]
			Ped<-Trim_ped[,(array_ped[s]): (array_ped[s]+((SNPcount*2)-1))]
			Ped_ready<-cbind(as.data.frame(SAMPLES),as.data.frame(Ped))
			
			#Create a list of two digits to be added, so the chromosomal segments can be sorted correctly
			segment_number<-sprintf("%02d",1:length(array_map))
			OUTPUT_map<-paste(WORKDIR,"/Chr",i,"_", segment_number[s],".map",sep="")
			OUTPUT_ped<-paste(WORKDIR,"/Chr",i,"_",segment_number[s],".ped",sep="")
			write.table(Map,OUTPUT_map,quote=F,row.names=F,col.names=F,sep="\t")
			write.table(Ped_ready,OUTPUT_ped,quote=F,row.names=F,col.names=F,sep="\t")
			
			#Make a list of the Starting and End SNP for each chromosomal segment of XSNPs
			Segments_positions<-Chr_map[c(array_map[s]),]
			Segments_info<-cbind(paste("Chr",i,"_",segment_number[s],sep=""),c("start"), Chr_map[array_map[s],])
			Segments_info2<-cbind(paste("Chr",i,"_",segment_number[s],sep=""),c("end"), Chr_map[c(array_map[s]+(SNPcount -1)),])
			SEGMENTS_INFORMATION<-cbind(Segments_info, Segments_info2)
			
			write.table(SEGMENTS_INFORMATION,paste(WORKDIR,"/Segments_info.txt",sep=""),append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")
			
			write.table(Map,paste(WORKDIR,"/map_cm_",SNPcount,"snpALL.map",sep=""),append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")
			
			}
		
}

#Run plink using the command line
#run Plink_run.sh . Create an "output" directory
#in the terminal make a list of all the file names without extension 
#find . -type f -name "*ped" | sed 's/\.\///g' | sed 's/.ped//g' >./List_gral.txt
