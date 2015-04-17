#This code calculate Fst comparing all populations to each other in a pairwise manner

rm(list=ls())
library(gtools)
library(hierfstat)

INPUT_file<-read.table("/home/morrellp/gonzales/Mohsen/New_Analysis_04162015/NorthAm_3453_2367snp_geneticorder.txt",header=T)

INPUT_infile[1:10,1:10]

Samples_info<-INPUT_infile[,c(1:4)]

GENOTYPES_POLY1<-INPUT_infile[,-c(1:4)]

#Convert AB genotypes to AA=1, BB=2,AB=NA

ConvertABto12<-function(dat){
	dat[dat == 'AA'] <-"1"
	dat[dat == 'BB'] <-"2"
	dat[dat =='AB'] <-"NA"
	
	return(dat)
}

Genotypes_numeric2<-as.data.frame(apply(GENOTYPES_POLY1,2, ConvertABto12))


Genotypes_numeric2[1:10,1:10]
#Get SNPs names
SNPs<-colnames(Genotypes_numeric2)

#Put together samples info and Genotypes 

INPUT<-cbind(as.data.frame(Samples_info),as.data.frame(Genotypes_numeric2))
#Count how many samples are in each breeding program
INPUT[1:10,1:10]
table(INPUT$BP_Code)

GENOTYPES_POLY<-INPUT[,-c(1:4)]
#Select the position of each program by row-type in the INPUT file

####
#AC  BA  MN  MT  N2  N6  OR  UT  VT  WA 
#382 293 384 315 384 384 298 377 327 381

####

#Select the position of each program by row-type in the INPUT file
AB2<-which(INPUT[,1] == 'AB' & INPUT[,2]== '2')
AB6<-which(INPUT[,1] == 'AB' & INPUT[,2]== '6')

BA2<-which(INPUT[,1] == 'BA' & INPUT[,2]== '2')
BA6<-which(INPUT[,1] == 'BA' & INPUT[,2]== '6')

MN6<-which(INPUT[,1] == 'MN' & INPUT[,2]== '6')

MT2<-which(INPUT[,1] == 'MT' & INPUT[,2]== '2')

N2<-which(INPUT[,1] == 'N2' & INPUT[,2]== '2')
N6<-which(INPUT[,1] == 'N6' & INPUT[,2]== '6')

OR6<-which(INPUT[,1] == 'OR' & INPUT[,2]== '6')

UT2<-which(INPUT[,1] == 'UT' & INPUT[,2]== '2')
UT6<-which(INPUT[,1] == 'UT' & INPUT[,2]== '6')

VT6<-which(INPUT[,1] == 'VT' & INPUT[,2]== '6')

WA2<-which(INPUT[,1] == 'WA' & INPUT[,2]== '2')
WA6<-which(INPUT[,1] == 'WA' & INPUT[,2]== '6')

##### 1. Comparing all pops to all pops


POPULATIONS<-c("AB2","AB6","BA2","BA6","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")

Fst_matrix<-matrix(NA,nrow=15,ncol=15)
colnames(Fst_matrix)<-c("AB2","AB6","BA2","BA6","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")
rownames(Fst_matrix)<-c("AB2","AB6","BA2","BA6","MN6","MT2","N2","N6","OR2","OR6","UT2","UT6","VT6","WA2","WA6")
	
	### Pop1
	focal_pop<-get(POPULATIONS[1])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 1:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[1,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
    
    ####Pop2
    focal_pop<-get(POPULATIONS[2])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 2:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[2,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
     
     #####pop3
     focal_pop<-get(POPULATIONS[3])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 3:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[3,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
      
      ####pop4
      focal_pop<-get(POPULATIONS[4])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 4:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[4,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
	
	###pop5
	focal_pop<-get(POPULATIONS[5])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 5:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[5,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
      ###pop6
	focal_pop<-get(POPULATIONS[6])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 6:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[6,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }

 ###pop7
	focal_pop<-get(POPULATIONS[7])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 7:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[7,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
         ###pop8
	focal_pop<-get(POPULATIONS[8])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 8:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[8,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
         ###pop9
	focal_pop<-get(POPULATIONS[9])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 9:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[9,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
       
        ###pop10
	focal_pop<-get(POPULATIONS[10])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 10:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[10,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
       
        ###pop11
	focal_pop<-get(POPULATIONS[11])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 11:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[11,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
       
        ###pop12
	focal_pop<-get(POPULATIONS[12])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 12:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[12,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
         ###pop13
	focal_pop<-get(POPULATIONS[13])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 13:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[13,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
         ###pop14
	focal_pop<-get(POPULATIONS[14])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 14:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[14,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
        ###pop15
	focal_pop<-get(POPULATIONS[15])
	FOCAL<-GENOTYPES_POLY[c(focal_pop),]
	#Each pop to all others
	for(o in 15:length(POPULATIONS)){
		
	other_pop<-get(POPULATIONS[o])
	TOTAL<-GENOTYPES_POLY[c(other_pop),]
	
	#Get all genotypes together, first the sub-population (0) then the total poupulation (1)
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))

	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	head(levels)
	#look at the Fst at each locus

     Fst_matrix[15,o]<-varcomp.glob(levels,loci,diploid=FALSE)$F[1,1]  
        }
        
      write.table(Fst_matrix,"/home/morrellp/gonzales/Mohsen/New_Analysis_04162015/Fst_matrix_ALL_pairwise.txt",quote=F,row.names=T,col.names=T,sep="\t")

