#title           :Folded_SFS_SinglePop.R  
#description     :Create a folded SFS.
#author		 :T. Kono
#========================================================================

#   What is the missing data value?
missing <- "NN"

#   A function to return whether or not a call is heterozygous
is_het <- function(x) {
    alleles <- unlist(strsplit(x, ""))
    counts <- unique(alleles)
    if( length(counts) != 1 ) {
        return(TRUE)
    }
    else {
        return(FALSE)
    }
}

#   A function to calculate the minor allele frequency
maf <- function(x) {
    #   remove missing genotype calls
    genotype <- x[x != missing]
    #   and heterozygous calls
    hets <- sapply(genotype, is_het)
    genotype <- genotype[!hets]
    #   Then count up each class
    g_counts <- table(genotype)
    #   If the table has no data, then we ditch it and return NA
    if(length(g_counts) == 0) {
        return(NA)
    }
    #   The names of the table are the different genotypes
    states <- names(g_counts)
    #   What is the frequency of A and B?
    #   First count up the AA, AB and BB
    AA <- g_counts[states[1]]
    BB <- g_counts[states[2]]
    #   This is a bit clumsy, but if any genotype is not present, it will take
    #   a value of NA. We set it to 0 in this case
    if(is.na(AA)) {
        AA <- 0
    }
    if(is.na(BB)) {
        BB <- 0
    }
    p <- AA / sum(g_counts)
    q <- BB / sum(g_counts)
    #   If either are 1, then return NA, as the marker is monomorphic
    #   in this population
    if(p == 1 || q == 1) {
        return(NA)
    }
    else {
        #   Return the minimum of those
        return(min(p, q))
    }
}

#   A function to create bins by MAF and count how many markers fall in
#   those bins
create_bins <- function(mafs, classes=10) {
    #   Create the breakpoints of the bins. We add 1 to the number of desired
    #   classes, since the number of _spaces_ between breakpoints is one less
    #   than the number of breakpoints.
    freq_classes <- seq(0, 0.5, length.out=classes + 1)
    #   Then put the MAFs into these classes
    binned_mafs <- cut(mafs, breaks=freq_classes, include.lowest=TRUE)
    #   Then convert it all to proportion
    binned_mafs_prop <- table(binned_mafs) / length(mafs)
    return(binned_mafs_prop)
}

#   Take arguments
args <- commandArgs(TRUE)
#   Which file is which?
geno <- args[1]

#   Read in the data
geno_dat <- read.table(geno, header=TRUE, row.names=1)

#   Calculate the SFS data for set 1
geno_mafs <- apply(geno_dat, 2, maf)
geno_sfs <- create_bins(geno_mafs)

#   Open a PDF file for output
filename <- gsub(".txt", "_SFS.pdf", args[1])
#   Get the breeding program and row type
partition <- unlist(strsplit(args[1], "[.]"))[1]
pdf(file=filename, 8, 6)
#   Plot the data
plt <- barplot(
    t(geno_sfs),
    ylim=c(0, 0.5),
    beside=TRUE,
    axisnames=F,
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main=paste("Folded SFS for", partition),
    col="black")
#   Create the labels
labels <- row.names(geno_sfs)
#   The plt object contains a matrix of positions of the bars
#   we will average along columns to get the positions
at <- apply(plt, 2, mean)
#   Stick them at the right place
axis(
    side=1,
    at=at,
    labels=labels,
    font=1,
    cex.axis=0.75)
dev.off()
