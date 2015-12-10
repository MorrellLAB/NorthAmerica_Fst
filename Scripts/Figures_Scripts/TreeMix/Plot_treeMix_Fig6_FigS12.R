#Plot treeMix results
#output names without suffix

#Figure 6 and S12.

#Source code from the treeMix package
source("plotting_funcs.R")


poporder<-c("Landraces",	"Idaho(2)",	"BuschAg(2)",	"BuschAg(Int.)",	"Montana(2)",	"NorthDakota(2)",	"Utah(2)",	"Washington(2)",	"Idaho(6)"	,"BuschAg(6)",	"Minnesota(6)"	,"NorthDakota(6)"	,"Utah(6)",	"Washington(6)",	"Oregon(2)"	,"Oregon(6)"	,"Virginia(6)")
#poporder<-c("AB2s", "BA2s" ,"BAI2s" ,"MT2s" ,"N2s", "UT2s", "WA2s" ,"AB6s", "BA6s", "MN6s" ,"N6s" ,"UT6s" ,"WA6s" ,"OR2w" ,"OR6w" ,"VT6w")
write.table(poporder,"poporder",quote=F,row.names=F,col.names=F,sep="\t")


#Draw the tree with arrows indicating the populations involved in introgression
pdf("~/mig0.pdf",width=5,height=5)
par(xpd=F)
plot_tree("out_run")
dev.off()

#Plot residuals
pdf("~/mig0r.pdf",width=4.2,height=4.2)
par(mar=c(7,7,2,2))
plot_resid("out_run", "poporder",cex=0.9)
dev.off()