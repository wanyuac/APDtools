# This is a subordinate function of SPDcomp.sh. It works in the similar way as
# merge_dist_tab.R, but input tables of reference distances and query distances
# are all produced by Bandage. Output format of this function is the same as
# merge_dist_tab.R.
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
#   Rscript --vanilla merge_dr_ds.R --drs 'dr1.tsv,dr2.tsv' --ds 'd1.tsv,d2.tsv'
# --out "strain1__Dm.tsv"
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 4 Nov 2018; the latest edition: 5 Nov 2018

# Read arguments ###############
library(optparse)

options = list(
    # input files
    make_option("--dr", type = "character", action = "store", default = NULL,
                help = "Comma-delimted paths to compiled reference distances"),
    make_option("--d", type = "character", action = "store", default = NULL,
                help = "Comma-delimited paths to files of compiled distances produced by compile_dists.py"),

    # filters for distances
    make_option("--d_max", type = "integer", action = "store", default = 300e3,
                help = "Maximum shortest-path distance (bp) to be filtered for. Set d_max < 0 to turn the filter off. [default = %default]"),
    make_option("--n_max", type = "integer", action = "store", default = 10,
                help = "Maximum node number in distance paths to be filtered for. Set n_max <= 0 to turn this filter off. [default = %default]"),

    # output arguments
    make_option("--out", type = "character", action = "store", default = "./merged_APDs.tsv",
                help = "Output directory [default = %default]")
)

opt.parser = OptionParser(option_list = options)
args = parse_args(opt.parser)


# Function definitions ###############
.reorderQueryNames <- function(d) {
    # A subordinate function of importTables. It sort query1 and query2 to simplify
    # merge and deduplication of data frames.
    for (i in 1 : nrow(d)) {
        qs <- d[i, c("query1", "query2")]
        d[i, c("query1", "query2")] <- qs[order(qs, decreasing = FALSE)]
    }

    return(d)
}

.filterDistances <- function(D, d_max, n_max) {
    # Filter the measurements where specified for reducing the time of merging
    if (d_max >= 0) {
        D <- subset(D, distance <= d_max)
        n <- nrow(D)
        if (n > 0) {
            print(paste(n, "rows passed the filter for a maximum distance of",
                        d_max, "bp.", sep = " "))
        } else {
            stop(paste0("Error: no distance measurements left after filtering for d_max (",
                        d_max, " bp)."))
        }
    }

    if (n_max > 0) {
        D <- subset(D, node_number <= n_max)
        n <- nrow(D)
        if (n > 0) {
            print(paste0(n, " rows passed the filter for a maximum node number of ",
                         n_max, "."))
        } else {
            stop(paste0("Error: no distance measurements left after filtering for n_max (",
                        n_max, ")."))
        }
    }

    return(D)
}

importTables <- function(tabs, d_max = -1, n_max = -1) {
    ds <- data.frame(query1 = character(0), query2 = character(0), distance = integer(0),
                     node_number = integer(0), stringsAsFactors = FALSE)
    for (t in tabs) {
        d <- read.delim(file = t, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        d <- d[, c("query1", "query2", "distance", "node_number")]
        d <- .filterDistances(D = d, d_max = d_max, n_max = n_max)
        ds <- rbind.data.frame(ds, .reorderQueryNames(d), stringsAsFactors = FALSE)
    }
    ds <- ds[!duplicated(ds[, c("query1", "query2")]), ]  # de-duplication

    return(ds)
}


# Main program ###############
# Setting up
ds <- strsplit(x = args$d, split = ",", fixed = TRUE)[[1]]
drs <- strsplit(x = args$dr, split = ",", fixed = TRUE)[[1]]
print(paste(length(drs), "tables of reference distances.", sep = " "))
print(paste(length(ds), "tables of query distances.", sep = " "))

# Import data
Dr <- importTables(tabs = drs)
Dr <- Dr[, c("query1", "query2", "distance")]
names(Dr)[3] <- "dr"

D <- importTables(tabs = ds, d_max = args$d_max, n_max = args$n_max)
names(D)[3] <- "d"

# Merge
Dm <- merge(x = Dr, y = D, by = c("query1", "query2"), all = FALSE, sort = FALSE)  # query1, query2, dr, d, node_number
Dm <- Dm[, c("query1", "query2", "d", "node_number", "dr")]
Dm$error <- Dm$d - Dm$dr

# Save result
write.table(Dm, file = args$out, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
print("Done.")
