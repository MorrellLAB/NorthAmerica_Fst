#   Create a derived site frequency spectrum from a file containing a single
#   column of derived allele frequencies.

#   These are the "Set1" Colors from the RColorBrewer package with 8 levels
colors <- c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#FFFF33",
    "#A65628",
    "#F781BF"
    )

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

#   Read in our matrices
ab_2row <- read.table("AB_2Row.SFS.txt", header=FALSE)
ab_2row <- as.numeric(ab_2row$V1)
ab_6row <- read.table("AB_6Row.SFS.txt", header=FALSE)
ab_6row <- as.numeric(ab_6row$V1)
bai_2row <- read.table("BAI_2Row.SFS.txt", header=FALSE)
bai_2row <- as.numeric(bai_2row$V1)
ba_2row <- read.table("BA_2Row.SFS.txt", header=FALSE)
ba_2row <- as.numeric(ba_2row$V1)
ba_6row <- read.table("BA_6Row.SFS.txt", header=FALSE)
ba_6row <- as.numeric(ba_6row$V1)
mn_6row <- read.table("MN_6Row.SFS.txt", header=FALSE)
mn_6row <- as.numeric(mn_6row$V1)
mt_2row <- read.table("MT_2Row.SFS.txt", header=FALSE)
mt_2row <- as.numeric(mt_2row$V1)
nd_2row <- read.table("N2_2Row.SFS.txt", header=FALSE)
nd_2row <- as.numeric(nd_2row$V1)
nd_6row <- read.table("N6_6Row.SFS.txt", header=FALSE)
nd_6row <- as.numeric(nd_6row$V1)
or_2row <- read.table("OR_2Row.SFS.txt", header=FALSE)
or_2row <- as.numeric(or_2row$V1)
or_6row <- read.table("OR_6Row.SFS.txt", header=FALSE)
or_6row <- as.numeric(or_6row$V1)
ut_2row <- read.table("UT_2Row.SFS.txt", header=FALSE)
ut_2row <- as.numeric(ut_2row$V1)
ut_6row <- read.table("UT_6Row.SFS.txt", header=FALSE)
ut_6row <- as.numeric(ut_6row$V1)
vt_6row <- read.table("VT_6Row.SFS.txt", header=FALSE)
vt_6row <- as.numeric(vt_6row$V1)
wa_2row <- read.table("WA_2Row.SFS.txt", header=FALSE)
wa_2row <- as.numeric(wa_2row$V1)
wa_6row <- read.table("WA_6Row.SFS.txt", header=FALSE)
wa_6row <- as.numeric(wa_6row$V1)
all_data <- read.table("All_SFS.txt", header=FALSE)
all_data <- as.numeric(all_data$V1)

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

#   Put the 2-row breeding programs together
sfs_data_2row <- cbind(
    ab_2row,
    bai_2row,
    ba_2row,
    mt_2row,
    nd_2row,
    or_2row,
    ut_2row,
    wa_2row
    )
#   These are the labels for the 2-row partitions
labels_2row <- c(
    "Idaho (2)",
    "Busch Ag. (Int.)",
    "Busch Ag. (2)",
    "Montana (2)",
    "North Dakota (2)",
    "Oregon (2)",
    "Utah (2)",
    "Washington (2)"
    )
filename_2row <- "Unfolded_SFS_2Row.pdf"
pdf(file=filename_2row, 10, 6)
#   Plot the data
plt <- barplot(
    t(sfs_data_2row),
    ylim=c(0, 0.70),
    beside=TRUE,
    axisnames=F,
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    col=colors,
    space=c(0, 2)
    )
#   Get labels for the x-axis tick marks
labels <- row.names(sfs_data_2row)
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
#   Add a legend
legend(
    "topright",
    labels_2row,
    fill=colors,
    cex=1.0
    )
dev.off()
#   Plot the same data, but do it in black and white. This may be a bit intense
#   pattern-wise, but it should work...
filename_2row <- "Unfolded_SFS_2Row_BW.pdf"
pdf(file=filename_2row, 10, 6)
#   Plot the data
plt <- barplot(
    t(sfs_data_2row),
    ylim=c(0, 0.70),
    beside=TRUE,
    axisnames=F,
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    col=c("#FFFFFF", "#BBBBBB", "#444444", "#000000"),
    space=c(0, 2)
    )
barplot(
    t(sfs_data_2row),
    ylim=c(0, 0.70),
    beside=TRUE,
    axisnames=F,
    density=c(8, 8, 8, 8, 32, 32, 32, 32),
    angle=c(0, 45, 0, 135),
    space=c(0, 2),
    col=c("#333333", "#333333", "#CCCCCC", "#CCCCCC"),
    add=TRUE
    )
#   Get labels for the x-axis tick marks
labels <- row.names(sfs_data_2row)
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
#   Add a legend
legend(
    "topright",
    labels_2row,
    fill=c("#FFFFFF", "#BBBBBB", "#444444", "#000000"),
    cex=1.25,
    ncol=2
    )
legend(
    "topright",
    labels_2row,
    density=c(8, 8, 8, 8, 32, 32, 32, 32),
    angle=c(0, 45, 0, 135),
    fill=c("#333333", "#333333", "#CCCCCC", "#CCCCCC"),
    cex=1.25,
    ncol=2,
    )
dev.off()


#   And the same with the 6-row programs
sfs_data_6row <- cbind(
    ab_6row,
    ba_6row,
    mn_6row,
    nd_6row,
    or_6row,
    ut_6row,
    vt_6row,
    wa_6row
    )
labels_6row <- c(
    "Idaho (6)",
    "Busch Ag. (6)",
    "Minnesota (6)",
    "North Dakota (6)",
    "Oregon (6)",
    "Utah (6)",
    "Virgina (6)",
    "Washington (6)"
    )
filename_6row <- "Unfolded_SFS_6row.pdf"
pdf(file=filename_6row, 10, 6)
#   Plot the data
plt <- barplot(
    t(sfs_data_6row),
    ylim=c(0, 0.70),
    beside=TRUE,
    axisnames=F,
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    col=colors,
    space=c(0, 2)
    )
#   Get labels for the x-axis tick marks
labels <- row.names(sfs_data_6row)
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
#   Add a legend
legend(
    "topright",
    labels_6row,
    fill=colors,
    cex=1.0
    )
dev.off()
#   Plot the same data, but do it in black and white. This may be a bit intense
#   pattern-wise, but it should work...
filename_6row <- "Unfolded_SFS_6Row_BW.pdf"
pdf(file=filename_6row, 10, 6)
#   Plot the data
plt <- barplot(
    t(sfs_data_6row),
    ylim=c(0, 0.70),
    beside=TRUE,
    axisnames=F,
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    col=c("#FFFFFF", "#BBBBBB", "#444444", "#000000"),
    space=c(0, 2)
    )
barplot(
    t(sfs_data_6row),
    ylim=c(0, 0.70),
    beside=TRUE,
    axisnames=F,
    density=c(8, 8, 8, 8, 32, 32, 32, 32),
    angle=c(0, 45, 0, 135),
    space=c(0, 2),
    col=c("#333333", "#333333", "#CCCCCC", "#CCCCCC"),
    add=TRUE
    )
#   Get labels for the x-axis tick marks
labels <- row.names(sfs_data_6row)
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
#   Add a legend
legend(
    "topright",
    labels_6row,
    fill=c("#FFFFFF", "#BBBBBB", "#444444", "#000000"),
    cex=1.25,
    ncol=2
    )
legend(
    "topright",
    labels_6row,
    density=c(8, 8, 8, 8, 32, 32, 32, 32),
    angle=c(0, 45, 0, 135),
    fill=c("#333333", "#333333", "#CCCCCC", "#CCCCCC"),
    cex=1.25,
    ncol=2,
    )
dev.off()

###
#   This block is for creating a HUGE derived SFS for each partition together
###
#   And put them together
# sfs_data <- cbind(
#     ab_2row,
#     ab_6row,
#     bai_2row,
#     ba_2row,
#     ba_6row,
#     mn_6row,
#     mt_2row,
#     nd_2row,
#     nd_6row,
#     or_2row,
#     or_6row,
#     ut_2row,
#     ut_6row,
#     vt_6row,
#     wa_2row,
#     wa_6row,
#     all_data
#     )
# filename <- "Unfolded_SFS_AllBP.pdf"
# pdf(file=filename, 10, 6)
# #   Plot the data
# plt <- barplot(
#     t(sfs_data),
#     ylim=c(0, 0.70),
#     beside=TRUE,
#     axisnames=F,
#     xlab="Derived Allele Frequency",
#     ylab="Proportion",
#     main="Derived Site Frequency Spectrum",
#     col=colors,
#     space=c(0, 2)
#     )
# #   Create the labels
# labels <- row.names(sfs_data)
# #   The plt object contains a matrix of positions of the bars
# #   we will average along columns to get the positions
# at <- apply(plt, 2, mean)
# #   Stick them at the right place
# axis(
#     side=1,
#     at=at,
#     labels=labels,
#     font=1,
#     cex.axis=0.75)
# #   Add a legend
# legend(
#     "topright",
#     c(
#         "AB 2-Row",
#         "AB 6-Row",
#         "BAI",
#         "BA 2-Row",
#         "BA 6-Row",
#         "MN 6-Row",
#         "MT 2-Row",
#         "ND 2-Row",
#         "ND 6-Row",
#         "OR 2-Row",
#         "OR 6-Row",
#         "UT 2-Row",
#         "UT 6-Row",
#         "VT 6-Row",
#         "WA 2-Row",
#         "WA 6-Row",
#         "All Samples"
#         ),
#     fill=colors,
#     cex=0.6
#     )
# dev.off()
