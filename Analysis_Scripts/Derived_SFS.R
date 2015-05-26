#   Create a derived site frequency spectrum from a file containing a single
#   column of derived allele frequencies.


#   A function to create bins by MAF and count how many markers fall in
#   those bins
create_bins <- function(mafs, classes=10) {
    #   Create the breakpoints of the bins. We add 1 to the number of desired
    #   classes, since the number of _spaces_ between breakpoints is one less
    #   than the number of breakpoints.
    freq_classes <- seq(0, 1, length.out=classes + 1)
    #   Then put the MAFs into these classes
    binned_mafs <- cut(mafs, breaks=freq_classes, include.lowest=TRUE)
    #   Then convert it all to proportion
    binned_mafs_prop <- table(binned_mafs) / length(mafs)
    return(binned_mafs_prop)
}

#   Take arguments
args <- commandArgs(TRUE)
daf <- read.table(args[1], header=FALSE)
daf <- as.numeric(daf$V1)

#   Read in the data
sfs <- create_bins(daf)

filename <- gsub(".txt", "_SFS.pdf", args[1])
#   Get the breeding program and row type
partition <- unlist(strsplit(args[1], "[.]"))[1]
pdf(file=filename, 8, 6)
#   Plot the data
plt <- barplot(
    t(sfs),
    ylim=c(0, max(sfs)),
    beside=TRUE,
    axisnames=F,
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main=paste("Unfolded SFS for", partition),
    col="black")
#   Create the labels
labels <- row.names(sfs)
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
