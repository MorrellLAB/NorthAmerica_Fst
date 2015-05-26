#   Create a folded SFS of two different filtering criteria to explor why the
#   numbers for average pairwise diversity change so much.

#   What is the missing data value?
missing <- "NA"

#   A function to calculate the minor allele frequency
maf <- function(x) {
    #   remove missing genotype calls
    genotype <- x[x != missing]
    #   Then count up each class
    g_counts <- table(genotype)
    #   If the table has no data, then we ditch it and return NA
    if(length(g_counts) == 0) {
        return(NA)
    }
    #   What is the frequency of A and B?
    #   First count up the AA, AB and BB
    AA <- g_counts["AA"]
    AB <- g_counts["AB"]
    BB <- g_counts["BB"]
    #   This is a bit clumsy, but if any genotype is not present, it will take
    #   a value of NA. We set it to 0 in this case
    if(is.na(AA)) {
        AA <- 0
    }
    if(is.na(AB)) {
        AB <- 0
    }
    if(is.na(BB)) {
        BB <- 0
    }
    p <- (2 * AA + AB) / (2 * sum(g_counts))
    q <- (2 * BB + AB) / (2 * sum(g_counts))
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
set1 <- args[1]
set2 <- args[2]

#   Read in the data
set1_dat <- read.table(set1, header=T)
set2_dat <- read.table(set2, header=T)

#   Calculate the SFS data for set 1
set1_mafs <- apply(set1_dat[,6:ncol(set1_dat)], 2, maf)
set1_sfs <- create_bins(set1_mafs)
#   And for set 2
set2_mafs <- apply(set2_dat[,6:ncol(set2_dat)], 2, maf)
set2_sfs <- create_bins(set2_mafs)

#   Stick them together into a data frame for plotting
sfs_data <- as.data.frame(cbind(set1_sfs, set2_sfs))
#   Open a PDF file for output
filename <- gsub(".txt", "_SFS.pdf", args[2])
#   Get the breeding program and row type
partition <- unlist(strsplit(args[2], "[.]"))[1]
pdf(file=filename, 10, 6)
#   Plot the data
plt <- barplot(
    t(sfs_data),
    ylim=c(0, max(sfs_data)),
    beside=TRUE,
    axisnames=F,
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main=paste("Folded SFS for", partition),
    col=c("black", "grey50"))
#   Create the labels
labels <- row.names(sfs_data)
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
#   And create a legend
legend(
    inset=0,
    cex=1,
    "topright",
    c("May 14 Dataset", "Het. Filter, Corrected Names"),
    fill=c("black", "grey50"))
dev.off()
