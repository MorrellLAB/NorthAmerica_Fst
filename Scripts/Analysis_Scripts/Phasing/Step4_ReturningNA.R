#Title           :Step4_ReturningNA.R 
#Description     :By default fastPHASE imputes missing data. Here we convert imputed values 
#				  back to missing.
#Author		 	 :A. Poets 
#Note		     :Requires output from Step3_Processing_output.sh
#========================================================================================

#The input file is the product of Step3. 
rm(list=ls())

#For each chromosome phased, replace missingess for NA

for (c in 1:7){
#fastPhased output modified to match format in prephase files
PHASE_OUT <-read.table(paste("~/MYoutput_",c,".txt",sep=""),header=F,row.names=1)

INPUT<-read.table(paste("~/prephase_files/NorthAm_prephase_chr_",c,".txt",sep=""),header=T,row.names=1)


#Now go SNP by SNP, identify '?' and replace those in the OUTPUT with 'NA'

for (s in 1:(dim(PHASE_OUT)[2])){

	MISSING<-which(INPUT[,s] == '?')
	
	PHASE_OUT[c(MISSING),s] <-NA 
}


colnames(PHASE_OUT)<-colnames(INPUT)

write.table(PHASE_OUT,paste("~/NorthAm_phased_chr",c,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")

}