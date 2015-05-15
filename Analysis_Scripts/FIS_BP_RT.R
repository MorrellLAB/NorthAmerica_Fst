#   R script to calculate per-breeding program and per-row type inbreeding
#   coefficients for the Barley CAP FST project
#   Thomas Kono
#       Jan 23 2015
#       Saint Paul, MN

#   Note: Requires diploid calls

#   What is missing data called in our data frame?
missing <- "NA"

#   Take arguments
args <- commandArgs(TRUE)

#   Define some functions to help us calculate FIS for each group
CalcFIS <- function(genotype)
{
    #   Calculates inbreeding coefficients as FIS given a vector of diploid
    #   genotype calls.
    #   First we drop missing calls, since they should not contribute to our
    #   calculation
    genotype <- genotype[!genotype == missing]
    #   Since each genotype is called A/C, this gets easy!
    #   Count up the numbers in each class
    genotype_counts <- table(genotype)
    #   Do some simple checks here:
    #   if genotype is length 0, then return NA
    if(!length(genotype) > 1)
    {
        return(NA)
    }
    #   if there are no 'AB' calls, we have no hets, and FIS is 1
    if(!"AB" %in% genotype)
    {
        return(1)
    }
    #   We want the frequency of the alleles
    p <- (2*genotype_counts["AA"] + genotype_counts["AB"]) / (2*sum(genotype_counts))
    q <- (2*genotype_counts["BB"] + genotype_counts["AB"]) / (2*sum(genotype_counts))
    #   FIS is defined as 1 - [Ho/He]
    He <- 2*p*q
    Ho <- genotype_counts["AB"] / sum(genotype_counts)
    FIS <- 1 - (Ho/He)
    return(FIS)
}

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
all_breeding <- read.table(args[1], header=T, stringsAsFactors=F)

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


### The following statements will throw some warnings, but that is just because
#   not every breeding program has 2-row and 6-row lines
#   In the cases where there are no lines, the data frame will be empty, which 
#   results in NaN results for FIS.
### This is not unexpected, and does not ruin our analysis.
#   We get NAs in our FIS calculations when there are genotype classes that do
#   not exist when we expect them to. For instance, if there are heterozygotes
#   but one of the homozygote classes is 0
#   Slice off the first 5 columns, since they are not genotype data
AB.2row.FIS <- mean(apply(AB.2row[,6:ncol(AB.2row)], 2, function(x) CalcFIS(x)), na.rm=T)
AB.6row.FIS <- mean(apply(AB.6row[,6:ncol(AB.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
BA.2row.FIS <- mean(apply(BA.2row[,6:ncol(BA.2row)], 2, function(x) CalcFIS(x)), na.rm=T)
BA.6row.FIS <- mean(apply(BA.6row[,6:ncol(BA.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
BAI.FIS <- mean(apply(BAI[,6:ncol(BA.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
MN.6row.FIS <- mean(apply(MN.6row[,6:ncol(MN.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
MT.2row.FIS <- mean(apply(MT.2row[,6:ncol(MT.2row)], 2, function(x) CalcFIS(x)), na.rm=T)
N2.FIS <- mean(apply(N2[,6:ncol(N2)], 2, function(x) CalcFIS(x)), na.rm=T)
N6.FIS <- mean(apply(N6[,6:ncol(N6)], 2, function(x) CalcFIS(x)), na.rm=T)
OR.2row.FIS <- mean(apply(OR.2row[,6:ncol(OR.2row)], 2, function(x) CalcFIS(x)), na.rm=T)
OR.6row.FIS <- mean(apply(OR.6row[,6:ncol(OR.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
UT.2row.FIS <- mean(apply(UT.2row[,6:ncol(UT.2row)], 2, function(x) CalcFIS(x)), na.rm=T)
UT.6row.FIS <- mean(apply(UT.6row[,6:ncol(UT.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
VT.6row.FIS <- mean(apply(VT.6row[,6:ncol(VT.6row)], 2, function(x) CalcFIS(x)), na.rm=T)
WA.2row.FIS <- mean(apply(WA.2row[,6:ncol(WA.2row)], 2, function(x) CalcFIS(x)), na.rm=T)
WA.6row.FIS <- mean(apply(WA.6row[,6:ncol(WA.6row)], 2, function(x) CalcFIS(x)), na.rm=T)

#   This is just for the sake of printing it nicely
print(
    list(
        AB.2row=AB.2row.FIS,
        AB.6Row=AB.6row.FIS,
        BA.2row=BA.2row.FIS,
        BA.6Row=BA.6row.FIS,
        BAI=BAI.FIS,
        MN.6Row=MN.6row.FIS,
        MT.2row=MT.2row.FIS,
        N2=N2.FIS,
        N6=N6.FIS,
        OR.2row=OR.2row.FIS,
        OR.6Row=OR.6row.FIS,
        UT.2row=UT.2row.FIS,
        UT.6Row=UT.6row.FIS,
        VT.6Row=VT.6row.FIS,
        WA.2row=WA.2row.FIS,
        WA.6Row=WA.6row.FIS
        )
    )
