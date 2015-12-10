#   Create a folded SFS of two different filtering criteria to explor why the
#   numbers for average pairwise diversity change so much.

#   What is the missing data value?
missing <- "NN"

#   These are the "Set1" colors from RColorBrewer, slightly modified for
#   plotting the number of classes that we have to here.
colors <- c(
    "#E41A1C",  #   AB2
    "#FD4D4F",  #   AB6
    "#377EB8",  #   BAI
    "#67AEE8",  #   BA2
    "#97DEFD",  #   BA6
    "#4DAF4A",  #   MN6
    "#984EA3",  #   MT2
    "#FF7F00",  #   ND2
    "#FFBF55",  #   ND6
    "#DDDD33",  #   OR2
    "#FFFF99",  #   OR6
    "#A65628",  #   UT2
    "#D9995D",  #   UT6
    "#F781BF",  #   VT6
    "#666666",  #   WA2
    "#AAAAAA",  #   WA6
    "#000000"  #   All
    )

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

#   Read in our matrices
ab_2row <- read.table("AB_2Row.txt", header=TRUE, row.names=1)
ab_6row <- read.table("AB_6Row.txt", header=TRUE, row.names=1)
bai_2row <- read.table("BAI_2Row.txt", header=TRUE, row.names=1)
ba_2row <- read.table("BA_2Row.txt", header=TRUE, row.names=1)
ba_6row <- read.table("BA_6Row.txt", header=TRUE, row.names=1)
mn_6row <- read.table("MN_6Row.txt", header=TRUE, row.names=1)
mt_2row <- read.table("MT_2Row.txt", header=TRUE, row.names=1)
nd_2row <- read.table("N2_2Row.txt", header=TRUE, row.names=1)
nd_6row <- read.table("N6_6Row.txt", header=TRUE, row.names=1)
or_2row <- read.table("OR_2Row.txt", header=TRUE, row.names=1)
or_6row <- read.table("OR_6Row.txt", header=TRUE, row.names=1)
ut_2row <- read.table("UT_2Row.txt", header=TRUE, row.names=1)
ut_6row <- read.table("UT_6Row.txt", header=TRUE, row.names=1)
vt_6row <- read.table("VT_6Row.txt", header=TRUE, row.names=1)
wa_2row <- read.table("WA_2Row.txt", header=TRUE, row.names=1)
wa_6row <- read.table("WA_6Row.txt", header=TRUE, row.names=1)
all_data <- read.table("All_Samples.txt", header=TRUE, row.names=1)

#   Calculate the minor allele frequencies
ab_2row <- apply(ab_2row, 2, maf)
ab_6row <- apply(ab_6row, 2, maf)
bai_2row <- apply(bai_2row, 2, maf)
ba_2row <- apply(ba_2row, 2, maf)
ba_6row <- apply(ba_6row, 2, maf)
mn_6row <- apply(mn_6row, 2, maf)
mt_2row <- apply(mt_2row, 2, maf)
nd_2row <- apply(nd_2row, 2, maf)
nd_6row <- apply(nd_6row, 2, maf)
or_2row <- apply(or_2row, 2, maf)
or_6row <- apply(or_6row, 2, maf)
ut_2row <- apply(ut_2row, 2, maf)
ut_6row <- apply(ut_6row, 2, maf)
vt_6row <- apply(vt_6row, 2, maf)
wa_2row <- apply(wa_2row, 2, maf)
wa_6row <- apply(wa_6row, 2, maf)
all_data <- apply(all_data, 2, maf)

#   Next, create the bins
ab_2row <- create_bins(ab_2row)
ab_6row <- create_bins(ab_6row)
bai_2row <- create_bins(bai_2row)
ba_2row <- create_bins(ba_2row)
ba_6row <- create_bins(ba_6row)
mn_6row <- create_bins(mn_6row)
mt_2row <- create_bins(mt_2row)
nd_2row <- create_bins(nd_2row)
nd_6row <- create_bins(nd_6row)
or_2row <- create_bins(or_2row)
or_6row <- create_bins(or_6row)
ut_2row <- create_bins(ut_2row)
ut_6row <- create_bins(ut_6row)
vt_6row <- create_bins(vt_6row)
wa_2row <- create_bins(wa_2row)
wa_6row <- create_bins(wa_6row)
all_data <- create_bins(all_data)

#   And put them together
sfs_data <- cbind(
    ab_2row,
    ab_6row,
    bai_2row,
    ba_2row,
    ba_6row,
    mn_6row,
    mt_2row,
    nd_2row,
    nd_6row,
    or_2row,
    or_6row,
    ut_2row,
    ut_6row,
    vt_6row,
    wa_2row,
    wa_6row,
    all_data
    )

filename <- "Folded_SFS_AllBP.pdf"
pdf(file=filename, 10, 6)
#   Plot the data
plt <- barplot(
    t(sfs_data),
    ylim=c(0, 0.5),
    beside=TRUE,
    axisnames=F,
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main="Folded Site Frequency Spectrum",
    col=colors)
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
    "topright",
    c(
        "AB 2-Row",
        "AB 6-Row",
        "BAI",
        "BA 2-Row",
        "BA 6-Row",
        "MN 6-Row",
        "MT 2-Row",
        "ND 2-Row",
        "ND 6-Row",
        "OR 2-Row",
        "OR 6-Row",
        "UT 2-Row",
        "UT 6-Row",
        "VT 6-Row",
        "WA 2-Row",
        "WA 6-Row",
        "All Samples"
        ),
    fill=colors,
    cex=0.6
    )
dev.off()
