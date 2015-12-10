#   R script to separate the large genotyping matrix into individual files,
#   based on breeding program and row type
#       Thomas Kono
#       April 16, 2015
#       Saint Paul, MN

#   Note: Requires diploid calls

#   What is missing data called in our data frame?
missing <- "NA"

#   Take arguments
args <- commandArgs(TRUE)

#   Rows are individuals and columns (from 6 until the end) are markers
#   First Eight columns are
#       Breeding Program
#           AB BA BAI MN MT N2 N6 OR UT VT WA
#       Sample ID
#       D1 (??)
#       D2 (??)
#       D3 (??)
#       Row Type (2 or 6)
#       CAP Year
#       Growth Habit (Winter or Spring)
all_breeding <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE)

#   Separate them by breeding program and row type
AB.2row <- all_breeding[which(all_breeding$Program == "AB" & all_breeding$Row_type == 2),]
AB.6row <- all_breeding[which(all_breeding$Program == "AB" & all_breeding$Row_type == 6),]

BA.2row <- all_breeding[which(all_breeding$Program == "BA" & all_breeding$Row_type == 2),]
BA.6row <- all_breeding[which(all_breeding$Program == "BA" & all_breeding$Row_type == 6),]
BAI <- all_breeding[which(all_breeding$Program == "BAI"),]

MN.2row <- all_breeding[which(all_breeding$Program == "MN" & all_breeding$Row_type == 2),]
MN.6row <- all_breeding[which(all_breeding$Program == "MN" & all_breeding$Row_type == 6),]

MT.2row <- all_breeding[which(all_breeding$Program == "MT" & all_breeding$Row_type == 2),]
MT.6row <- all_breeding[which(all_breeding$Program == "MT" & all_breeding$Row_type == 6),]

#   Already separated by RT
N2 <- all_breeding[which(all_breeding$Program == "N2"),]
N6 <- all_breeding[which(all_breeding$Program == "N6"),]

OR.2row <- all_breeding[which(all_breeding$Program == "OR" & all_breeding$Row_type == 2),]
OR.6row <- all_breeding[which(all_breeding$Program == "OR" & all_breeding$Row_type == 6),]

UT.2row <- all_breeding[which(all_breeding$Program == "UT" & all_breeding$Row_type == 2),]
UT.6row <- all_breeding[which(all_breeding$Program == "UT" & all_breeding$Row_type == 6),]

VT.2row <- all_breeding[which(all_breeding$Program == "VT" & all_breeding$Row_type == 2),]
VT.6row <- all_breeding[which(all_breeding$Program == "VT" & all_breeding$Row_type == 6),]

WA.2row <- all_breeding[which(all_breeding$Program == "WA" & all_breeding$Row_type == 2),]
WA.6row <- all_breeding[which(all_breeding$Program == "WA" & all_breeding$Row_type == 6),]

#   Write the tables out
write.table(AB.2row, "AB_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(AB.6row, "AB_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(BA.2row, "BA_2Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(BA.6row, "BA_6Row.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(BAI, "BAI.txt", quote=FALSE, sep="\t", row.names=FALSE)
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
