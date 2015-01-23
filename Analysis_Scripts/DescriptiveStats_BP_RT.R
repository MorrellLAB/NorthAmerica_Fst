#   R script to calculate per-breeding program and per-row type diversity stats
#   for the Barley CAP FST project
#   Thomas Kono
#       Jan 23 2015
#       Saint Paul, MN

#   Define some functions that will help us calculate pairwise similarity
#       This one calculates allele counts
AlleleCounts <- function(genotype)
{
    #   We remove missing calls, since they shouldn't contribute to pairwise
    #   similarity
    genotype <- genotype[!is.na(genotype)]
    #   We get the allele counts
    allele_counts <- table(genotype)
    #   Return a list with the major and minor allele counts as elements
    return(list(major=max(allele_counts), minor=min(allele_counts)))
}
#       This one calculates pairwise similarity given allele counts
PairwiseSimilarity <- function(major_count, minor_count)
{
    #   The way to calculate pairwise similarity for a single marker is
    #   [(NumAllele1 choose 2) + (NumAllele2 choose 2)]/(AllAlleles choose 2)
    #   Essentially, the number of similar comparisons (allele1 to allele1
    #   and allele2 to allele2) divided by the number of comparisons total
    Allele1 <- choose(major_count, 2)
    Allele2 <- choose(minor_count, 2)
    AllAlleles <- choose(major_count + minor_count, 2)
    pairwise <- (Allele1 + Allele2)/AllAlleles
    return(pairwise)
}

#   Rows are individuals and columns (from 6 until the end) are markers
#   First five columns are
#       Breeding Program
#           AC BA MN MT N2 N6 OR UT VT WA
#       Sample ID
#       Row Type (2 or 6)
#       CAP Year
#       Growth Habit (Winter or Spring)
all_breeding <- read.table('BreedingGenotypes_homo_QC.txt', header=T)

#   Separate them by breeding program and row type
AC.2row <- all_breeding[which(all_breeding$BP_Code == "AC" & all_breeding$Row_Type == 2),]
AC.6row <- all_breeding[which(all_breeding$BP_Code == "AC" & all_breeding$Row_Type == 6),]

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


### The following statements will throw some warnings, but that is just because
#   not every breeding program has 2-row and 6-row lines
#   In the cases where there are no lines, the data frame will be empty, which 
#   results in NaN results for pairwise similarity.
### This is not unexpected, and does not ruin our analysis.
AC.2row.PS <- mean(apply(AC.2row[,6:ncol(AC.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
AC.6row.PS <- mean(apply(AC.6row[,6:ncol(AC.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
BA.2row.PS <- mean(apply(BA.2row[,6:ncol(BA.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
BA.6row.PS <- mean(apply(BA.6row[,6:ncol(BA.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
MN.2row.PS <- mean(apply(MN.2row[,6:ncol(MN.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
MN.6row.PS <- mean(apply(MN.6row[,6:ncol(MN.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
MT.2row.PS <- mean(apply(MT.2row[,6:ncol(MT.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
MT.6row.PS <- mean(apply(MT.6row[,6:ncol(MT.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
N2.PS <- mean(apply(N2[,6:ncol(N2)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
N6.PS <- mean(apply(N6[,6:ncol(N6)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
OR.2row.PS <- mean(apply(OR.2row[,6:ncol(OR.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
OR.6row.PS <- mean(apply(OR.6row[,6:ncol(OR.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
UT.2row.PS <- mean(apply(UT.2row[,6:ncol(UT.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
UT.6row.PS <- mean(apply(UT.6row[,6:ncol(UT.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
VT.2row.PS <- mean(apply(VT.2row[,6:ncol(VT.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
VT.6row.PS <- mean(apply(VT.6row[,6:ncol(VT.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
WA.2row.PS <- mean(apply(WA.2row[,6:ncol(WA.2row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
WA.6row.PS <- mean(apply(WA.6row[,6:ncol(WA.6row)], 2, function(x) PairwiseSimilarity(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)

#   This is just for the sake of printing it nicely
print(
    list(
        AC.2row=AC.2row.PS,
        AC.6Row=AC.6row.PS,
        BA.2row=BA.2row.PS,
        BA.6Row=BA.6row.PS,
        MN.2row=MN.2row.PS,
        MN.6Row=MN.6row.PS,
        MT.2row=MT.2row.PS,
        MT.6Row=MT.6row.PS,
        N2=N2.PS,
        N6=N6.PS,
        OR.2row=OR.2row.PS,
        OR.6Row=OR.6row.PS,
        UT.2row=UT.2row.PS,
        UT.6Row=UT.6row.PS,
        VT.2row=VT.2row.PS,
        VT.6Row=VT.6row.PS,
        WA.2row=WA.2row.PS,
        WA.6Row=WA.6row.PS
        )
    )
