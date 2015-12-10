###Estimate and plot (Figure S4) derived SFS and joint two SFS comparing spike types or growth habit types.

rm(list=ls())

###INPUT FILES: Can be found in the Dataset folder.
### load marker file (row-samples; column-markers). Use diploid nucleotide calls. 
dat<-read.table("Barley_NorthAm_QC_ACTG_no_duplicates_or.txt",header=T,row.names=1)

#     X12_30969 X12_10420 X11_10895 X11_11223 X11_21354
#AB_1        CC        AA        GG        GG        CC
#AB_2        CC        AA        GG        GG        CC
#AB_3        CC        AA        GG        GG        CC


sam<-read.table("samples_information.txt",header=T)
### load sample file (row type and growth habit). File with four columns: Alias (or name), AccessionID, Program, Row-type, Growth_habit
#	Alias       Name Program Row_type Growth_habit
#	AB_1 04AB029-22      AB        2            S
#	AB_2 04AB003-59      AB        2            S
#	AB_3 04AB098-51      AB        2            S

### load ancestral state for each marker. This is table S3 from Fang et al., 2014.
ancestral<-read.table("Fang_et_al_TableS3.txt",header=T)

#     SNP Morex_State Bulbosum_States
#11_10005           A               G
#12_30341           G               G
#12_30718           G              CG


##FILES FORMATING:
#make ancestral states diploid 
ancestral_bul<-as.character(ancestral[,3])
ancestral_bul[as.character(ancestral_bul)=="A"]<-"AA"
ancestral_bul[as.character(ancestral_bul)=="C"]<-"CC"
ancestral_bul[as.character(ancestral_bul)=="T"]<-"TT"
ancestral_bul[as.character(ancestral_bul)=="G"]<-"GG"

#Mark as missing data those SNPs that are trans-specific between barley and H. bulbosum.
ancestral_bul[(ancestral_bul!="AA") & (ancestral_bul!="CC") & (ancestral_bul!="TT") & (ancestral_bul!="GG")] <- NA

ancestral_dip<-as.data.frame(cbind(as.character(ancestral[,1]), ancestral_bul))
anc<-ancestral_dip[order(ancestral_dip[,1]),]

#Turn genotypes to have samples in columns and SNPs in rows
dat2 <- t(dat)
rownames(dat2)<-substr(rownames(dat2), 2,9)
dat3 <- dat2[rownames(dat2) %in% anc[,1], ]
anc2 <- anc[anc[,1] %in% rownames(dat2), ]
dat4 <- dat3[order(rownames(dat3)), ]
dat5 <- t(dat4)

### convert the ancestral state to 0, derived state to 2, heterozygotes to NA
for (i in 1:dim(anc2)[1]) {
  dat5[,i][which(dat5[,i] == as.character(anc2[i,2]))] <- "0"
  }

### the ancestral state is neither the two alleles on some of the SNPs, remove these SNPs
noanc <- NULL
for (i in 1:dim(anc2)[1]) {
  noanc[i] <- "0" %in% dat5[,i]
  }

dat6 <- dat5[,-which(noanc=="FALSE")]

dat6[dat6=="AA"]<-"2"
dat6[dat6=="TT"]<-"2"
dat6[dat6=="CC"]<-"2"
dat6[dat6=="GG"]<-"2"
dat6[dat6=="NN"]<-"NA"
dat6[(dat6!="0") & (dat6!="2") & (dat6!="NA")] <- "1"

dat7<-cbind(rownames(dat6), dat6)
colnames(dat7)[1] <- "Alias"
dat8 <- merge(sam[,-2], dat7, by="Alias")

tworow <- subset(dat8, dat8[,3]=="2")
sixrow <- subset(dat8, dat8[,3]=="6")
spring <- subset(dat8, dat8[,4]=="S")
winter <- subset(dat8, dat8[,4]=="W")

tworow2 <- tworow[,c(-1,-2,-3,-4)]
sixrow2 <- sixrow[,c(-1,-2,-3,-4)]
spring2 <- spring[,c(-1,-2,-3,-4)]
winter2 <- winter[,c(-1,-2,-3,-4)]

## SFS CALCULATION:
### calculate derived allele frequency
n0e <- apply(tworow2==0,2,sum,na.rm=T)
n1e <- apply(tworow2==1,2,sum,na.rm=T)
n2e <- apply(tworow2==2,2,sum,na.rm=T)
ne <- n0e + n1e + n2e
maf2 <- ((2*n2e)+n1e)/(2*ne)

n0w <- apply(sixrow2==0,2,sum,na.rm=T)
n1w <- apply(sixrow2==1,2,sum,na.rm=T)
n2w <- apply(sixrow2==2,2,sum,na.rm=T)
nw <- n0w + n1w + n2w
maf6 <- ((2*n2w)+n1w)/(2*nw)

n0e <- apply(spring2==0,2,sum,na.rm=T)
n1e <- apply(spring2==1,2,sum,na.rm=T)
n2e <- apply(spring2==2,2,sum,na.rm=T)
ne <- n0e + n1e + n2e
mafspring <- ((2*n2e)+n1e)/(2*ne)

n0w <- apply(winter2==0,2,sum,na.rm=T)
n1w <- apply(winter2==1,2,sum,na.rm=T)
n2w <- apply(winter2==2,2,sum,na.rm=T)
nw <- n0w + n1w + n2w
mafwinter <- ((2*n2w)+n1w)/(2*nw)

##JOINT SFS PLOTTING:
### joint SFS heatmap for two and six spike type
vec1<-maf2
vec2<-maf6

space=1
total=length(vec1)

xs=hist(vec1,plot=F)$breaks
bklongx=space*length(xs)
xs=hist(vec1,plot=F,breaks=bklongx)$breaks

ys=hist(vec2,plot=F)$breaks
bklongy=space*length(ys)
ys=hist(vec2,plot=F,breaks=bklongy)$breaks
comb=data.frame(vec1,vec2)
library(grDevices)

c=matrix(0,(length(xs)-1),(length(ys)-1))
for( i in 1:(length(xs)-1)){
for( j in 1:(length(ys)-1)){

ay=subset(comb,comb$vec1>=xs[i])
bee=subset(ay,ay$vec1<xs[i+1])
cee=subset(bee,bee$vec2>=ys[j])
d=subset(cee,cee$vec2<ys[j+1])
c[i,j]=length(d$vec2)/total
}
}

### filled.contour(c, color = terrain.colors, plot.axes={})
pdf("Two-six_jSFS.pdf",width=6,height=5)
library(gplots)
rgb.palette<-colorRampPalette(c("white","red","black"),space="rgb")
filled.contour(c, col = rgb.palette(30), xlab="Two-row barley", ylab="Six-row barley")
dev.off()

#Plot Joint SFS for spring and winter types.
vec2<-mafspring
vec1<-mafwinter

space=1
total=length(vec1)

xs=hist(vec1,plot=F)$breaks
bklongx=space*length(xs)
xs=hist(vec1,plot=F,breaks=bklongx)$breaks

ys=hist(vec2,plot=F)$breaks
bklongy=space*length(ys)
ys=hist(vec2,plot=F,breaks=bklongy)$breaks
comb=data.frame(vec1,vec2)
library(grDevices)

c=matrix(0,(length(xs)-1),(length(ys)-1))
for( i in 1:(length(xs)-1)){
for( j in 1:(length(ys)-1)){

ay=subset(comb,comb$vec1>=xs[i])
bee=subset(ay,ay$vec1<xs[i+1])
cee=subset(bee,bee$vec2>=ys[j])
d=subset(cee,cee$vec2<ys[j+1])
c[i,j]=length(d$vec2)/total
}
}

### filled.contour(c, color = terrain.colors, plot.axes={})
pdf("Winter-spring_jSFS.pdf",width=6,height=5)
library(gplots)
rgb.palette<-colorRampPalette(c("white","red","black"),space="rgb")
filled.contour(c, col = rgb.palette(30), xlab="Winter barley", ylab="Spring barley")
dev.off()
