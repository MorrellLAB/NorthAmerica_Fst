

rm<-list(ls())

##For Spring 6 row

#Call mean Fst values for all 100000 iterations run over each SNP
cvalue<-read.table("NorthAm_c_value_global_spring6.txt",header=T)

#Plot density curve of c-values (each c-value is the result of one iteration across all SNPs)

pdf("Spring_6_Nicholson.pdf")
plot(density(cvalue[,1]),col="blue",ylim=c(0,150),xlim=c(0,0.4),lwd=2, xlab="c-value",main="Marginal posterior densities for Spring 6-rows")
lines(density(cvalue[,2]),col="dark green",lwd=2)
lines(density(cvalue[,3]),col="maroon",lwd=2)
lines(density(cvalue[,4]),col="purple",lwd=2)
lines(density(cvalue[,5]),col="cyan",lwd=2)
lines(density(cvalue[,6]),col="orange",lwd=2)
		
legend("topright",cex=0.9,legend=c("Idaho","Busch Agriculture", "Minnesota", "North Dakota" , "Utah", "Washington"),col=c("blue","dark green","maroon","purple","cyan","orange"),lty=c(rep(1,6)),lwd=2,title=expression(bold("Breeding Program")),title.adj=0.70)
dev.off()


##For Spring 2-row
rm(list=ls())
#Call mean Fst values for all 100000 iterations run over each SNP
cvalue<-read.table("NorthAm_c_value_global_spring2.txt",header=T)
head(cvalue)
#Plot density curve of c-values (each c-value is the result of one iteration across all SNPs)

pdf("Spring_2_Nicholson.pdf")
plot(density(cvalue[,1]),col="blue",ylim=c(0,200),xlim=c(0,0.4),lwd=2, xlab="c-value",main="Marginal posterior densities for Spring 2-rows")
lines(density(cvalue[,2]),col="dark green",lwd=2)
lines(density(cvalue[,3]),col="green",lwd=2)
lines(density(cvalue[,4]),col="gray",lwd=2)
lines(density(cvalue[,5]),col="purple",lwd=2)
lines(density(cvalue[,6]),col="cyan",lwd=2)
lines(density(cvalue[,7]),col="orange",lwd=2)

		
legend("topright",cex=0.9,legend=c("Idaho","Busch Ag.","Busch Ag.(Int.)", "Montana", "North Dakota" , "Utah", "Washington"),col=c("blue","dark green","GREEN","gray","purple","cyan","orange"),lty=c(rep(1,6)),lwd=2,title=expression(bold("Breeding Program")),title.adj=0.70)
dev.off()


##TABLE S10:
#Make a table of mean and std for all the partitions
sp6_mean<-read.table("NorthAm_muc_value_global_spring6.txt",header=T,row.names=1)
sp2_mean<-read.table("NorthAm_muc_value_global_spring2.txt",header=T,row.names=1)

sp6_std<-read.table("NorthAm_sdc_value_global_spring6.txt",header=T,row.names=1)
sp2_std<-read.table("NorthAm_sdc_value_global_spring2.txt",header=T,row.names=1)


MEAN<-rbind(as.data.frame(sp6_mean),as.data.frame(sp2_mean))
STD<-rbind(sp6_std,sp2_std)

TABLE_ready<-format(cbind(MEAN,STD),scientific=F)
colnames(TABLE_ready)<-c("Mean c-value","Standard Devidation")

row.names(TABLE_ready)<-c("University of Idaho in Aberdeen (6-row)","Bush Agr. Resources (6-row)","Univesity of Minnesota","North Dakota State University (6-row)","Utah State University (6-row)","Washington State University (6-row)","University of Idaho in Aberdeen (2-row)","Bush Agr. Resources (2-row)","Bush Agr. Resources (International)","Montana State University", "North Dakota State University (2-row)","Utah State University (2-row)","Washington State University (2-row)")

write.table(TABLE_ready,"~/Dropbox/Mohsen_Fst/Tables/Fst_Nicholson_mean_std.xls",quote=F,row.names=T,col.names=T,sep="\t")



