#Code to make histograms of Fst distributions.
#Input files are the resutls from Fst_six_two_row.R, Fst_winter_spring.r, Fst_BreedingPopulation.R

rm(list=ls())

#For GH and RT and BP FST
FST1<-read.table("Fst_winter_spring.txt",header=F)
FST2<-read.table("Fst_RT.txt",header=F)
FST3<-read.table("Fst_BP.txt",header=F)

MAIN<-c("Growth Habit","Row Type", "Breeding Program")
FILE<-c("Fst_GH","Fst_RT","Fst_BP")

for (i in 1:3){
	hfst<-get(paste("FST",i,sep=""))

	#get quantile at 97.5%
	extreme<-quantile(hfst[,2],0.975, na.rm=T)

	h<-hist(hfst[,2])
	clr <- ifelse(h$breaks < extreme, "grey", "red")[-length(h$breaks)]
	AVERAGE<-summary((hfst[,2]))[4]
	pdf(paste(FILE[i],"_dist.pdf",sep=""),width=7,height=5)
	plot(h,col=clr,main=MAIN[i], xlab = "Fst")
	legend("topright",c("97.5th percentile cutoff"),pch=15,col="red")
	abline(v=AVERAGE, col="black",lty=2)
	dev.off()
}