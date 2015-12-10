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
breeding <- read.table(args[1], header=T, stringsAsFactors=F)

### The following statements will throw some warnings, but that is just because
#   not every breeding program has 2-row and 6-row lines
#   In the cases where there are no lines, the data frame will be empty, which 
#   results in NaN results for FIS.
### This is not unexpected, and does not ruin our analysis.
#   We get NAs in our FIS calculations when there are genotype classes that do
#   not exist when we expect them to. For instance, if there are heterozygotes
#   but one of the homozygote classes is 0
#   Slice off the first 5 columns, since they are not genotype data
FIS <- mean(apply(breeding[,2:ncol(breeding)], 2, function(x) CalcFIS(x)), na.rm=T)
print(c(args[1], FIS))
