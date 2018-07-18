#!/usr/bin/env Rscript
# This script calculates reference distances based on feature tables of the same bacterial strain. It assums every genome under
# analysis is circular and it ignores pseudo genes.
#
# Usage:
#   Rscript calc_dr.R --strain [strain name] --features [feature_table.txt] --contig_sizes [a TSV file for genome sizes] --nofuzzy_start
# The name of the output directory must be the same as the strain name under the current working directory.
# This script does not work for any genomes having only a single feature, such as one CDS.
#
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
# First edition: 28/6/2016; the latest edition: 12 Dec 2017

library(optparse)

options = list(
    make_option(c("-s", "--strain"), type = "character", default = NULL, dest = "strain",
                help = "Strain name", metavar="character"),
    make_option(c("-f", "--features"), type = "character", default = "feature_table.tsv", dest = "features",
                help = "A tab-delimited table of genomic features [default= %default]", metavar = "character"),
    make_option(c("-c", "--contig_sizes"), type = "character", default = "contig_sizes.tsv", dest = "contig_sizes",
                help = "A two-column tab-delimited table (without header) for contig sizes (bp) [default= %default]",
                metavar = "character"),
    make_option(c("-n", "--nofuzzy_start"), type = "logical", action = "store_true", dest = "nofuzzy_start",
                help = "Flag it when the \"feature.location.nofuzzy_start\" attribute of BioPython is used.",
                metavar = "logical")
)

COMPLETENESS_STATUS <- c(complete = 0, incomplete.start = 1, incomplete.end = 2, incomplete.both = 3)

# Define a function that calculates reference distances for every molecule ####################
getCircularPhysicalDistance <- function(ft, genome_size, filename) {
    contig_name <- ft$Contig_name[1]
    n <- nrow(ft)
    
    # fill other cells
    for (i in 1 : (n - 1)) {
        for (j in (i + 1) : n) {
            row_i <- ft[i, ]  # extract a single row to increase the speed of locating variables
            row_j <- ft[j, ]
            fi <- row_i$Locus_tag
            fj <- row_j$Locus_tag
            s <- c(row_i$Start, row_j$Start)
            e <- c(row_i$End, row_j$End)
            len <- c(row_i$Length, row_j$Length)  # gene lengths
            R <- max(e) - min(s) + 1  # the range covered by two features clockwisely
            Dc <- R - sum(len)  # clockwise distance
            Da <- genome_size - R  # anticlockwise distance >= 0
            D <- min(Dc, Da)
            if (D < 0) {
                D <- 0  # assign a distance of zero for overlapping genes
            }
            write(paste(c(contig_name, fi, fj, D), collapse = "\t"), file = filename, append = TRUE)
        }
    }
}

# Read arguments ##################
opt_parser <- OptionParser(option_list = options)
args <- parse_args(opt_parser)
strain <- args$strain  # The first element must be the strain name.

ft <- read.delim(args$features, stringsAsFactors = FALSE)  # read the feature table

# remove pseudo genes
if ("Pseudo" %in% names(ft)) {
    print("Excluding pseudo genes from the feature table.")
    
    # Empty values of the pseudo column may be "" when there are "*" values or NA when all genes are not pseudo.
    ft$Pseudo[which(is.na(ft$Pseudo))] <- ""  # substitute NA's with empty characters
    ft <- subset(ft, Pseudo != "*")
    if (nrow(ft) == 0) {
        stop("Warning: this script does not work when all features are pseudo genes.\n")
    }
    ft <- ft[, -which(names(ft) == "Pseudo")]
}

# determine the sequence completeness of each feature
feature_sum <- nrow(ft)  # the number of features
ft <- cbind(ft, Completeness = integer(feature_sum))
for (j in 1 : feature_sum) {
    start <- ft$Start[j]
    end <- ft$End[j]
    is_start_incomplete <- substr(start, 1, 1) == "<"  # whether the 5'end is truncated
    is_end_incomplete <- substr(end, 1, 1) == ">"  # whether the 3'end is truncated
    if (is_start_incomplete) {
        ft$Completeness[j] <- ft$Completeness[j] + COMPLETENESS_STATUS[["incomplete.start"]]
        ft$Start[j] <- substr(start, 2, nchar(start))  # remove the ">" character from the beginning
    }
    if (is_end_incomplete) {
        ft$Completeness[j] <- ft$Completeness[j] + COMPLETENESS_STATUS[["incomplete.end"]]
        ft$End[j] <- substr(end, 2, nchar(end))  # remove the "<" sign from the beginning
    }
}
ft$Start <- as.integer(ft$Start)
ft$End <- as.integer(ft$End)

# deal with the nofuzzy_start option
if (args$nofuzzy_start) {
    print("Adding one to the start coordinates to correct the shift due to the nofuzzy_start option in BioPython.")
    ft$Start <- ft$Start + 1
}

# generate some summary statistics for the remaining features
contigs <- unique(ft$Contig_name)
feature_count <- table(ft$Contig_name)  # number of features per contig
contig_num <- length(contigs)
cat("There are", contig_num, "contigs and", feature_sum, "non-pseudo features for the strain", paste0(strain, ".\n"))
print("Saving the final feature table.")
write.table(ft, file = paste0(strain, "__nopseudo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# import contig sizes
cs_df <- read.delim(args$contig_sizes, header = FALSE, stringsAsFactors = FALSE)[, 1 : 2]
names(cs_df) <- c("Contig_name", "Length")
cs <- cs_df$Length
names(cs) <- cs_df$Contig_name  # a named integer vector
rm(cs_df)

# Extract features of each contig and calculate reference distances ##################
# initialise the output file
filename <- paste0(strain, "__dr.tsv")
write(paste(c("Contig_name", "Query1", "Query2", "Distance"), collapse = "\t"), file = filename)  # A newline character is attached automatically.

# append lines to the file
for (ctg in contigs) {
    if (feature_count[[ctg]] > 1) {  # A distance is possible when there are at least two loci in the same contig.
        if ((ctg %in% names(cs))) {
            print(paste0("Calculating reference distances in the contig ", ctg, "."))
            getCircularPhysicalDistance(ft = subset(ft, Contig_name == ctg),
                                        genome_size = cs[[ctg]],
                                        filename = filename)
        } else {
            stop(paste("Error: size of the contig ", ctg, "is unknown.\n", sep = " "))
        }
    }
}

print("Reference distances are successfully calculated.")
