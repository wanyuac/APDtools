# Compare measured allelic physical distances to reference distances and
# calculate the accuracy accordingly under a given error tolerance.
#
# Parameters
#    d, dr: distances measured in de novo assemblies and complete genomes respectively.
#           Each input is a tab-delimited file produced using the script compile_dists.py.
#
# This script is designed to make it easier for users to obtain accuracy in their
# distance measurements in simulation-validaation studies. As this script works
# for scientific research and does not emphasis security, please do not challenge
# this program using dodgy parameters.
#
# Example command line:
#   Rscript --vanilla d_vs_dr.R --d d.tsv --dr dr.tsv --d_max 250000 --n_max 5 \\
#   --e_max 500 --outdir dist/comp --basename d_vs_dr
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 12 May 2018; the latest edition: 12 May 2018

# Read arguments ###############
library(optparse)

options = list(
    # input files
    make_option("--d", type = "character", action = "store", default = NULL,
                help = "Path to compiled distances measured in genome assemblies"),
    make_option("--dr", type = "character", action = "store", default = NULL,
                help = "Path to compiled reference distances"),

    # filters for distances
    make_option("--d_max", type = "integer", action = "store", default = 300e3,
                help = "Maximum shortest-path distance (bp) to be filtered for. Set d_max < 0 to turn the filter off. [default = %default]"),
    make_option("--n_max", type = "integer", action = "store", default = 10,
                help = "Maximum node number in distance paths to be filtered for. Set n_max <= 0 to turn this filter off. [default = %default]"),
    make_option("--e_max", type = "integer", action = "store", default = 1e3,
                help = "Tolerance to absolute values of errors (bp). [default = %default]"),

    # output arguments
    make_option("--outdir", type = "character", action = "store", default = ".",
                help = "Output directory [default = %default]"),
    make_option("--basename", type = "character", action = "store",
                default = "d_vs_dr", help = "Basename of output files [default = %default]")
)

opt.parser = OptionParser(option_list = options)
args = parse_args(opt.parser)

# Function definitions ###############
mergeDrD <- function(D, Dr) {
    # Merge two data frames into a single one.
    # For robustness, this function does not assume query1 and query2 follow the
    # same permutation in D and Dr, although this is the case when D and Dr are
    # obtained using Bandage with the same file of allele sequences (queries).
    samples <- unique(D$sample)
    Dm <- NULL  # initialisation

    for (s in samples) {
        Dm_s <- subset(D, sample == s)  # rows of D to be merged with dr
        n <- nrow(Dm_s)
        Dm_s <- cbind.data.frame(Dm_s, dr = integer(n), stringsAsFactors = FALSE)  # initialise the result data frame
        for (i in 1 : n) {
            r <- Dm_s[i, c("query1", "query2")]  # take a row from Dm
            q1 <- r[["query1"]]
            q2 <- r[["query2"]]
            Dr_s <- subset(Dr, sample == s)  # The sample must be present in Dr in normal scenarios.
            j <- which((Dr_s$query1 == q1 & Dr_s$query2 == q2) | (Dr_s$query1 == q2 & Dr_s$query2 == q1))
            nj <- length(j)
            if (nj == 1) {
                Dm_s$dr[i] <- Dr_s$dr[j]  # get the reference distance
            } else if (nj == 0) {
                # This may happen when distances are unmeasurable in complete genomes
                # but become measurable in de novo assemblies because of additional
                # connections.
                Dm_s$dr[i] <- NA
                print(paste0("Warning: D contains queries (", q1, ", ", q2, ") that are not present in Dr of the sample ",
                             s, "."))

            } else {  # nj > 1: an erroneous condition
                stop(paste0("Error: queries in D have multiple matches to Dr of the sample ",
                            s, "."))
            }
        }
        Dm <- rbind.data.frame(Dm, Dm_s, stringsAsFactors = FALSE)
    }


    # remove rows where dr is absent
    Dm <- subset(Dm, !is.na(dr))

    return(Dm)
}

evaluateAccuracy <- function(Dm, e_max) {
    # Initialisation
    print(paste("Calculating accuracy under an error tolerance of", e_max, "bp.",
                sep = " "))
    samples <- unique(Dm$sample)
    sample_n <- length(samples)
    ns <- 1 : max(Dm$node_number)
    acc_colnames <- c("Sample", paste0("N", ns))
    acc <- matrix(NA, nrow = sample_n, ncol = length(acc_colnames))
    colnames(acc) <- acc_colnames
    rownames(acc) <- NULL
    acc[, "Sample"] <- samples

    # Calculate the accuracy given each number of nodes
    for (n in ns) {
        Dm_n <- subset(Dm, node_number <= n)
        vn <- character(sample_n)  # values for the column Nn
        for (i in 1 : sample_n) {
            s <- samples[i]
            Dm_n_s <- subset(Dm_n, sample == s)
            d_num <- nrow(Dm_n_s)
            if (d_num > 0) {
                acc_count <- sum(Dm_n_s$e_abs <= e_max)
                acc_perc <- round(acc_count / d_num * 100, digits = 2)
                vn[i] <- paste0(acc_perc, "% (", acc_count, "/", d_num, ")")
            } else {
                vn[i] <- NA
            }
        }
        acc[, paste0("N", n)] <- vn
    }
    acc <- as.data.frame(acc, stringsAsFactors = FALSE)

    return(acc)
}

# Main program ###############
# Import distances measured in de novo genome assemblies
D <- read.delim(file = args$d, stringsAsFactors = FALSE)
D <- D[, c("query1", "query2", "sample", "distance", "node_number", "source")]
D <- D[order(D$sample, D$node_number, decreasing = FALSE), ]
names(D)[4] <- "d"
n <- nrow(D)
print(paste(n, "distance measurements are found."))
print("Before filtering for upper bounds of distances and node numbers:")
print(paste("  Maximum distance across samples:", max(D$d), "bp.", sep = " "))
print(paste0("  Maximum node number :", max(D$node_number), "."))

# Filter the measurements where specified
if (args$d_max >= 0) {
    D <- subset(D, d <= args$d_max)
    n <- nrow(D)
    if (n > 0) {
        print(paste(n, "rows passed the filter for a maximum distance of",
                     args$d_max, "bp.", sep = " "))
    } else {
        stop(paste0("Error: no distance measurements left after filtering for d_max (",
                    args$d_max, " bp)."))
    }
}

if (args$n_max > 0) {
    D <- subset(D, node_number <= args$n_max)
    n <- nrow(D)
    if (n > 0) {
        print(paste0(n, " rows passed the filter for a maximum node number of ",
                     args$n_max, "."))
    } else {
        stop(paste0("Error: no distance measurements left after filtering for n_max (",
                    args$n_max, ")."))
    }
}

# Import reference distances
Dr <- read.delim(file = args$dr, stringsAsFactors = FALSE)
Dr <- Dr[, c("query1", "query2", "sample", "distance")]
names(Dr)[4] <- "dr"

# Merge two tables into a single one
Dm <- mergeDrD(D, Dr)
Dm$e_abs <- abs(Dm$d - Dm$dr)
Dm <- Dm[, c("query1", "query2", "sample", "dr", "d", "e_abs", "node_number", "source")]

# Calculate the accuracy for every sample
acc <- evaluateAccuracy(Dm, args$e_max)

# Write the result into a tab-delimited file
write.table(acc, file = paste0(args$outdir, "/", args$basename, ".tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
write.table(Dm, file = paste0(args$outdir, "/", args$basename, "__Dm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

print("Done.")
