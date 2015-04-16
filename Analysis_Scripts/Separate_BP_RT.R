#   R script to separate the large genotyping matrix into individual files,
#   based on breeding program and row type
#       Thomas Kono
#       April 16, 2015
#       Saint Paul, MN

#   Note: Requires diploid calls

#   What is missing data called in our data frame?
missing <- "NA"


#   Rows are individuals and columns (from 6 until the end) are markers
#   First Eight columns are
#       Breeding Program
#           AB BA MN MT N2 N6 OR UT VT WA
#       Sample ID
#       D1 (??)
#       D2 (??)
#       D3 (??)
#       Row Type (2 or 6)
#       CAP Year
#       Growth Habit (Winter or Spring)
all_breeding <- read.table('NORTH_AMERICAN_QC_3453samples_2396snp.txt', header=T, stringsAsFactors=F)

#   Separate them by breeding program and row type
AB.2row <- all_breeding[which(all_breeding$BP_Code == "AB" & all_breeding$Row_Type == 2),]
AB.6row <- all_breeding[which(all_breeding$BP_Code == "AB" & all_breeding$Row_Type == 6),]

BA.2row <- all_breeding[which(all_breeding$BP_Code == "BA" & all_breeding$Row_Type == 2),]
BA.6row <- all_breeding[which(all_breeding$BP_Code == "BA" & all_breeding$Row_Type == 6),]

MN.2row <- all_breeding[which(all_breeding$BP_Code == "MN" & all_breeding$Row_Type == 2),]
MN.6row <- all_breeding[which(all_breeding$BP_Code == "MN" & all_breeding$Row_Type == 6),]

MT.2row <- all_breeding[which(all_breeding$BP_Code == "MT" & all_breeding$Row_Type == 2),]
MT.6row <- all_breeding[which(all_breeding$BP_Code == "MT" & all_breeding$Row_Type == 6),]

#   Already separated by RT
N2 <- all_breeding[which(all_breeding$BP_Code == "N2"),]
N6 <- all_breeding[which(all_breeding$BP_Code == "N6"),]

OR.2row <- all_breeding[which(all_breeding$BP_Code == "OR" & all_breeding$Row_Type == 2),]
OR.6row <- all_breeding[which(all_breeding$BP_Code == "OR" & all_breeding$Row_Type == 6),]

UT.2row <- all_breeding[which(all_breeding$BP_Code == "UT" & all_breeding$Row_Type == 2),]
UT.6row <- all_breeding[which(all_breeding$BP_Code == "UT" & all_breeding$Row_Type == 6),]

VT.2row <- all_breeding[which(all_breeding$BP_Code == "VT" & all_breeding$Row_Type == 2),]
VT.6row <- all_breeding[which(all_breeding$BP_Code == "VT" & all_breeding$Row_Type == 6),]

WA.2row <- all_breeding[which(all_breeding$BP_Code == "WA" & all_breeding$Row_Type == 2),]
WA.6row <- all_breeding[which(all_breeding$BP_Code == "WA" & all_breeding$Row_Type == 6),]

#   Write the tables out
write.table(AB.2row, "AB_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(AB.6row, "AB_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(BA.2row, "BA_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(BA.6row, "BA_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(MN.6row, "MN_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(MT.2row, "MT_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(N2, "ND_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(N6, "ND_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(OR.2row, "OR_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(OR.6row, "OR_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(UT.2row, "UT_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(UT.6row, "UT_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(VT.6row, "VT_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(WA.2row, "WA_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(WA.6row, "WA_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
