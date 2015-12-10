#Analysis after Running Step1_PED_data_generation.R and  Step2_Make_Nicholson_FST_Inputs.py  
 
# Get Fst genome-wide (global)

rm(list=ls())

#install packages from URL
#install.packages("http://cran.r-project.org/src/contrib/popgen_1.0-3.tar.gz", method="libcurl")
library(popgen)
#Set up working directories.
#INPUTDIR where the output of Step2 is.
INPUTDIR<-c("~/Input/")
OUTPUTDIR<-c("~/output/")

PARTITIONS<-c("spring2",'spring6',"winter6")
for (i in 1:length(PARTITIONS)){
	N<-read.table(paste(INPUTDIR,"NorthAm_",PARTITIONS[i],"_N.txt",sep=""))
	X<-read.table(paste(INPUTDIR,"NorthAm_",PARTITIONS[i],"_X.txt",sep=""))
	NUMA<-read.table(paste(INPUTDIR,"NorthAm_",PARTITIONS[i],"_NUMA.txt",sep=""))
	POPS<-read.table(paste(INPUTDIR,"Samples_pops_",PARTITIONS[i],".clst",sep=""))
	
	L <- dim(NUMA)[1] #SNP number
	P <- length(table(POPS[,2])) # Number of populations analyzed
	
##CALCULATING NICHOLSON c values genome-wide:

RESULTS_POPDIV2 <- popdiv(L, P, as.matrix(NUMA), as.matrix(N), as.matrix(X), burnin = 1000, iter = 10000, m = 10, csd = 0.01, outc = TRUE)

c_value<-RESULTS_POPDIV2$C.sample


if (PARTITIONS[i] == 'spring6'){
colnames(c_value)<-c("AB_6","BA_6","MN_6","N_6","UT_6","WA_6")
}

if (PARTITIONS[i] == 'spring2'){
colnames(c_value)<-c("AB_2","BA_2","BAI_2","MT_2","N_2","UT_2","WA_2")
}

if (PARTITIONS[i] == 'winter6'){
colnames(c_value)<-c("OR_6","VT_6")
}

#posterior mean and std of fst
muc_value <-RESULTS_POPDIV2$muc
muc_value <-cbind(colnames(c_value),muc_value)

sdc_value<-RESULTS_POPDIV2$sdc
sdc_value <-cbind(colnames(c_value), sdc_value)

#Results genome-wide per each population
write.table(c_value,paste(OUTPUTDIR,"NorthAm_c_value_global_",PARTITIONS[i],".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
write.table(muc_value,paste(OUTPUTDIR,"NorthAm_muc_value_global_",PARTITIONS[i],".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
write.table(sdc_value,paste(OUTPUTDIR,"NorthAm_sdc_value_global_",PARTITIONS[i],".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}


##These results are plotted using Plot_nicholson.R located in the Plot_Scripts directory.