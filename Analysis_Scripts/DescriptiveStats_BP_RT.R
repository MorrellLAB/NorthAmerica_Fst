#   R script to calculate per-breeding program and per-row type diversity stats
#   for the Barley CAP FST project
#   Thomas Kono
#       Jan 23 2015
#       Saint Paul, MN

#   Define some functions that will help us calculate pairwise Difference
#       This one calculates allele counts
AlleleCounts <- function(genotype)
{
    #   We remove missing calls, since they shouldn't contribute to pairwise
    #   Difference
    genotype <- genotype[!is.na(genotype)]
    #   We get the allele counts
    allele_counts <- table(genotype)
    #   Return a list with the major and minor allele counts as elements
    return(list(major=max(allele_counts), minor=min(allele_counts)))
}
#       This one calculates pairwise Difference given allele counts
PairwiseDifference <- function(major_count, minor_count)
{
    #   The way to calculate pairwise Difference for a single marker is
    #   1 - [(NumAllele1 choose 2) + (NumAllele2 choose 2)]/(AllAlleles choose 2)
    #   Essentially, 1 minus the number of similar comparisons (allele1 to allele1
    #   and allele2 to allele2) divided by the number of comparisons total
    Allele1 <- choose(major_count, 2)
    Allele2 <- choose(minor_count, 2)
    AllAlleles <- choose(major_count + minor_count, 2)
    pairwise <-  1 - ((Allele1 + Allele2)/AllAlleles)
    return(pairwise)
}

#   These are diagnostic functions
#       % Missingness
PropMissing <- function(genotype)
{
    total_calls <- length(genotype)
    #   how many missing
    num_missing <- sum(is.na(genotype))
    prop_missing <- num_missing/total_calls
    return(prop_missing)
}
#       Minor allele frequency
MAF <- function(genotype)
{
    #   First remove missing data
    genotype <- genotype[!is.na(genotype)]
    #   check if they are monomorphic
    if(length(unique(genotype)) == 1)
    {
        return(0)
    }
    #   count them up!
    allele_counts <- table(genotype)
    minor_allele <- min(allele_counts)
    #   and calculate the frequency
    maf <- minor_allele / length(genotype)
    return(maf)
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

MN.6row <- all_breeding[which(all_breeding$BP_Code == "MN" & all_breeding$Row_Type == 6),]

MT.2row <- all_breeding[which(all_breeding$BP_Code == "MT" & all_breeding$Row_Type == 2),]

#   Already separated by RT
N2 <- all_breeding[which(all_breeding$BP_Code == "N2"),]
N6 <- all_breeding[which(all_breeding$BP_Code == "N6"),]

OR.2row <- all_breeding[which(all_breeding$BP_Code == "OR" & all_breeding$Row_Type == 2),]
OR.6row <- all_breeding[which(all_breeding$BP_Code == "OR" & all_breeding$Row_Type == 6),]

UT.2row <- all_breeding[which(all_breeding$BP_Code == "UT" & all_breeding$Row_Type == 2),]
UT.6row <- all_breeding[which(all_breeding$BP_Code == "UT" & all_breeding$Row_Type == 6),]

VT.6row <- all_breeding[which(all_breeding$BP_Code == "VT" & all_breeding$Row_Type == 6),]

WA.2row <- all_breeding[which(all_breeding$BP_Code == "WA" & all_breeding$Row_Type == 2),]
WA.6row <- all_breeding[which(all_breeding$BP_Code == "WA" & all_breeding$Row_Type == 6),]

### The following statements will throw some warnings, but that is just because
#   not every breeding program has 2-row and 6-row lines
#   In the cases where there are no lines, the data frame will be empty, which 
#   results in NaN results for pairwise Difference.
### This is not unexpected, and does not ruin our analysis.
AC.2row.PD <- mean(apply(AC.2row[,6:ncol(AC.2row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
AC.6row.PD <- mean(apply(AC.6row[,6:ncol(AC.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
BA.2row.PD <- mean(apply(BA.2row[,6:ncol(BA.2row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
BA.6row.PD <- mean(apply(BA.6row[,6:ncol(BA.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
MN.6row.PD <- mean(apply(MN.6row[,6:ncol(MN.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
MT.2row.PD <- mean(apply(MT.2row[,6:ncol(MT.2row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
N2.PD <- mean(apply(N2[,6:ncol(N2)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
N6.PD <- mean(apply(N6[,6:ncol(N6)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
OR.2row.PD <- mean(apply(OR.2row[,6:ncol(OR.2row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
OR.6row.PD <- mean(apply(OR.6row[,6:ncol(OR.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
UT.2row.PD <- mean(apply(UT.2row[,6:ncol(UT.2row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
UT.6row.PD <- mean(apply(UT.6row[,6:ncol(UT.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
VT.6row.PD <- mean(apply(VT.6row[,6:ncol(VT.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
WA.2row.PD <- mean(apply(WA.2row[,6:ncol(WA.2row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)
WA.6row.PD <- mean(apply(WA.6row[,6:ncol(WA.6row)], 2, function(x) PairwiseDifference(AlleleCounts(x)$major, AlleleCounts(x)$minor)), na.rm=T)

#   Now we will calculate the proportion of missingness
AC.2row.missing <- mean(apply(AC.2row[,6:ncol(AC.2row)], 2, function(x) PropMissing(x)))
AC.6row.missing <- mean(apply(AC.6row[,6:ncol(AC.6row)], 2, function(x) PropMissing(x)))
BA.2row.missing <- mean(apply(BA.2row[,6:ncol(BA.2row)], 2, function(x) PropMissing(x)))
BA.6row.missing <- mean(apply(BA.6row[,6:ncol(BA.6row)], 2, function(x) PropMissing(x)))
MN.6row.missing <- mean(apply(MN.6row[,6:ncol(MN.6row)], 2, function(x) PropMissing(x)))
MT.2row.missing <- mean(apply(MT.2row[,6:ncol(MT.2row)], 2, function(x) PropMissing(x)))
N2.missing <- mean(apply(N2[,6:ncol(N2)], 2, function(x) PropMissing(x)))
N6.missing <- mean(apply(N6[,6:ncol(N6)], 2, function(x) PropMissing(x)))
OR.2row.missing <- mean(apply(OR.2row[,6:ncol(OR.2row)], 2, function(x) PropMissing(x)))
OR.6row.missing <- mean(apply(OR.6row[,6:ncol(OR.6row)], 2, function(x) PropMissing(x)))
UT.2row.missing <- mean(apply(UT.2row[,6:ncol(UT.2row)], 2, function(x) PropMissing(x)))
UT.6row.missing <- mean(apply(UT.6row[,6:ncol(UT.6row)], 2, function(x) PropMissing(x)))
VT.6row.missing <- mean(apply(VT.6row[,6:ncol(VT.6row)], 2, function(x) PropMissing(x)))
WA.2row.missing <- mean(apply(WA.2row[,6:ncol(WA.2row)], 2, function(x) PropMissing(x)))
WA.6row.missing <- mean(apply(WA.6row[,6:ncol(WA.6row)], 2, function(x) PropMissing(x)))

#   And the average minor allele frequency
AC.2row.maf <- mean(apply(AC.2row[,6:ncol(AC.2row)], 2, function(x) MAF(x)))
AC.6row.maf <- mean(apply(AC.6row[,6:ncol(AC.6row)], 2, function(x) MAF(x)))
BA.2row.maf <- mean(apply(BA.2row[,6:ncol(BA.2row)], 2, function(x) MAF(x)))
BA.6row.maf <- mean(apply(BA.6row[,6:ncol(BA.6row)], 2, function(x) MAF(x)))
MN.6row.maf <- mean(apply(MN.6row[,6:ncol(MN.6row)], 2, function(x) MAF(x)))
MT.2row.maf <- mean(apply(MT.2row[,6:ncol(MT.2row)], 2, function(x) MAF(x)))
N2.maf <- mean(apply(N2[,6:ncol(N2)], 2, function(x) MAF(x)))
N6.maf <- mean(apply(N6[,6:ncol(N6)], 2, function(x) MAF(x)))
OR.2row.maf <- mean(apply(OR.2row[,6:ncol(OR.2row)], 2, function(x) MAF(x)))
OR.6row.maf <- mean(apply(OR.6row[,6:ncol(OR.6row)], 2, function(x) MAF(x)))
UT.2row.maf <- mean(apply(UT.2row[,6:ncol(UT.2row)], 2, function(x) MAF(x)))
UT.6row.maf <- mean(apply(UT.6row[,6:ncol(UT.6row)], 2, function(x) MAF(x)))
VT.6row.maf <- mean(apply(VT.6row[,6:ncol(VT.6row)], 2, function(x) MAF(x)))
WA.2row.maf <- mean(apply(WA.2row[,6:ncol(WA.2row)], 2, function(x) MAF(x)))
WA.6row.maf <- mean(apply(WA.6row[,6:ncol(WA.6row)], 2, function(x) MAF(x)))


#   This is just for the sake of printing it nicely
print(
    list(
        AC=structure(c(AC.2row.PD, AC.6row.PD, AC.2row.missing, AC.6row.missing, AC.2row.maf, AC.6row.maf), names=c("2rowPD", "6rowPD", "2rowMissing", "6rowMissing", "2rowAvgMAF", "6rowAvgMAF")),
        BA=structure(c(BA.2row.PD, BA.6row.PD, BA.2row.missing, BA.6row.missing, BA.2row.maf, BA.6row.maf), names=c("2rowPD", "6rowPD", "2rowMissing", "6rowMissing", "2rowAvgMAF", "6rowAvgMAF")),
        MN=structure(c(MN.6row.PD, MN.6row.missing,MN.6row.maf), names=c("6rowPD", "6rowMissing", "6rowAvgMAF")),
        MT=structure(c(MT.2row.PD, MT.2row.missing, MT.2row.maf), names=c("2rowPD", "2rowMissing", "2rowAvgMAF")),
        N2=structure(c(N2.PD, N2.missing, N2.maf), names=c("2rowPD", "2rowMissing", "2rowAvgMAF")),
        N6=structure(c(N6.PD, N6.missing, N6.maf), names=c("6rowPD", "6rowMissing", "6rowAvgMAF")),
        OR=structure(c(OR.2row.PD, OR.6row.PD, OR.2row.missing, OR.6row.missing, OR.2row.maf, OR.6row.maf), names=c("2rowPD", "6rowPD", "2rowMissing", "6rowMissing", "2rowAvgMAF", "6rowAvgMAF")),
        UT=structure(c(UT.2row.PD, UT.6row.PD, UT.2row.missing, UT.6row.missing, UT.2row.maf, UT.6row.maf), names=c("2rowPD", "6rowPD", "2rowMissing", "6rowMissing", "2rowAvgMAF", "6rowAvgMAF")),
        VT=structure(c(VT.6row.PD, VT.6row.missing, VT.6row.maf), names=c("6rowPD", "6rowMissing", "6rowAvgMAF")),
        WA=structure(c(WA.2row.PD, WA.6row.PD, WA.2row.missing, WA.6row.missing, WA.2row.maf, WA.6row.maf), names=c("2rowPD", "6rowPD", "2rowMissing", "6rowMissing", "2rowAvgMAF", "6rowAvgMAF"))
        )
    )
