#   Script to calculate per-SNP pairwise diversity and and plot it in order
#   of a genetic map. libsequence-based tools (at least for the pairwise
#   diversity calculations) remove missing sites and then calculate on the
#   remaining bases. We will implement the same methodology here.

#   Take arguments
args <- commandArgs(TRUE)
#   The first argument will be the genetic map
genmap <- args[1]
#   The second argument will be the genotyping matrix
genotypes <- args[2]
#   Set the missing values. Here we use heterozygotes and missing genotypes
missing <- c(NA, "AB")
#   Set the colors for plotting. These are the "Set1" colors from RColorBrewer
#   with seven levels.
# colors <- c(
#     "#a6cee3",
#     "#1f78b4",
#     "#b2df8a",
#     "#33a02c",
#     "#fb9a99",
#     "#e31a1c",
#     "#fdbf6f"
#     )
#   These are the colors that Ana uses in the FST plots
colors <- c(
    "red",
    "purple",
    "green",
    "blue",
    "yellow",
    "dark green",
    "orange")
#   Light and dark grey, for printing
# colors <- c(
#     "grey70",
#     "black",
#     "grey70",
#     "black",
#     "grey70",
#     "black",
#     "grey70"
#     )
#   And the chromosomes we wish to plot
chrom <- c(
    "1H",
    "2H",
    "3H",
    "4H",
    "5H",
    "6H",
    "7H"
    )

#   Define a function to calculate the pairwise difference
pairwise_div <- function(x) {
    #   First, remove missing data
    #   is.element(a, b) will return a vector of boolean values of length
    #   length(a), indicating whether or not a[i] is in b.
    #   This removes heterozygous calls and missing calls from our genotypes.
    calls <- x[!is.element(x, missing)]
    #   Count up the genotypic classes
    counts <- table(calls)
    #   If the length of the counts is 1, then the SNP is monomorphic, and
    #   pairwise diversity is 0.
    if(length(counts) == 1) {
        pairwise_div <- 0
    }
    else {
        #   Then get the minimum and max for major and minor alleles
        major <- max(counts)
        minor <- min(counts)
        #   The way to calculate pairwise *similarity* for a single SNP is
        #   [(MinorCount choose 2) + (MajorCount choose 2)] / (N.Ind choose 2)
        #   Essentially, the number of similarities due to major/major
        #   comparisons and minor/minor comparisons, divided by the total number
        #   of comparisons.
        minorsim <- choose(minor, 2)
        majorsim <- choose(major, 2)
        nindsim <- choose(length(calls), 2)
        pairwise_sim <- (minorsim + majorsim) / nindsim
        #   Pairwise div. is 1-similarity
        pairwise_div <- 1 - pairwise_sim
    }
    return(pairwise_div)
}

#   Define a function to calculate sliding window means
sliding_window_mean <- function(values, coords, winsize=10, winstep=5) {
    #   Create a vector of window start positions
    points <- seq(from=1, to=length(values) - winsize, by=winstep)
    #   Start two new vectors to grow that contain average diversity and
    #   the midpoint of the window
    means <- c()
    midpoints <- c()
    for(i in 1:length(points)) {
        #   If we are well within the middle of the vector, we still have a
        #   "whole" window to average
        if( (i + winsize) <= length(values) ) {
            means[i] <- mean(values[points[i]:points[i] + winsize])
            midpoints[i] <- mean(coords[points[i]:points[i]+winsize])
        }
        else {
            #   Otehrwise, we take from the window start to the end of the
            #   vector of values. An incomplete window.
            means[i] <- mean(values[points[i]:length(values)])
            midpoints[i] <- mean(values[points[i]:length(values)])
        }
    }
    #   Then put the vectors together into a data frame for plotting
    sliding_windows <- data.frame(
        Diversity=means,
        Coords=midpoints)
    return(sliding_windows)
}

#   Read in the genetic map
genetic_map <- read.table(genmap, header=TRUE, stringsAsFactors=FALSE)
#   And the genotypes
genotyping_matrix <- read.table(genotypes, header=TRUE, row.names=1, stringsAsFactors=FALSE)

#   Get the per-SNP pairwise diversity numbers for each SNP
pairwise_diversity <- apply(genotyping_matrix, 2, pairwise_div)
#   Start a new plot
fname <- gsub(".txt", "_SlidingWindow_PairwiseDiv.pdf", args[2])
partition <- unlist(strsplit(args[2], "[.]"))[1]
pdf(file=fname, 10, 5)
plot(
    c(0, max(genetic_map$Cumulative)),
    c(0, 0.55),
    type="n",
    xlab="Linkage Group",
    ylab="Average Pairwise Diversity",
    main=paste("Average Pairwise Diversity in ", partition, sep=""),
    xaxt="n"
    )
#   For each chromosome, get the markers that are on it and calculate sliding
#   window means
chr_meds <- c()
for(i in 1:length(chrom)) {
    chromosome <- chrom[i]
    #   Which markers are on this chromosome?
    chr_markers <- genetic_map[genetic_map$chromosome == chromosome,]
    #   And what are the diversity values?
    chr_div <- pairwise_diversity[chr_markers$SNP]
    #   Then, sort them according to the order in the geneic map
    chr_div <- chr_div[match(chr_markers$SNP, names(chr_div))]
    #   Stick them into a data frame
    chr_div <- data.frame(
        Position=chr_markers$Cumulative,
        Marker=chr_markers$SNP,
        Diversity=as.numeric(chr_div)
        )
    #   Remove missing data
    chr_div <- chr_div[!is.na(chr_div$Diversity),]
    chr_sliding_windows <- sliding_window_mean(
        chr_div$Diversity,
        chr_div$Position
        )
    #   Then, we make the 0 values as really small, but not 0, so they will
    #   work with log
    chr_sliding_windows$Diversity[chr_sliding_windows$Diversity == 0] <- 1e-4
    #   Throw the points on the plot
    #   Same graphical plotting parameters that Ana uses for FST plot
    lines(
        chr_sliding_windows$Coords,
        chr_sliding_windows$Diversity,
        lwd=1,
        col=colors[i]
        )
    # points(
    #     chr_sliding_windows$Coords,
    #     chr_sliding_windows$Diversity,
    #     pch=20,
    #     cex=0.75,
    #     col=colors[i]
    #     )
    #   Save the midpoint of the cumulative value for labeling the plot
    #   later
    chr_meds[i] <- median(chr_sliding_windows$Coords)
}
#   Label the x-axis with the linkage groups
axis(
    side=1,
    at=chr_meds,
    labels=chrom,
    font=3
    )
dev.off()
