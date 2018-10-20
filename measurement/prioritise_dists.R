# With this script, a user can specify weights of distance measurements from different sources and filter out measurements of low weights if
# they also have high-weight ones.
# Usage example:
#   Rscript --vanilla prioritise_dists.R -d distances.tsv -w source_weights.tsv -o filtered_dists.tsv -c 16
# Column names of the weight table: source and weight. The larger the weight, the more superior the distance is. For example:
#   source  weight
#   finished_genome  3
#   contig    2
#   graph  1
# Author: Yu Wan (wanyuac@gmail.com)
# Development history: 23/11/2016, 23/5/2017
# Licence: GPL v2.1

library(optparse)
library(parallel)

# Obtain arguments from the command line ###############
options = list(
    make_option(c("-d", "--dists"), type = "character", default = NULL, help = "Input distance table in tab-delimited format", metavar="character"),
    make_option(c("-w", "--weights"), type = "character", default = "dist_weights.txt",
                help = "A two-column tab-delimited file about weights of sources: [default= %default]",
                metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = "filtered_dists.tsv",
                help = "Output file name [default= %default]", metavar="character"),
    make_option(c("-c", "--cores"), type = "integer", default = 0,
                help = "How many cores to be used [default= %default: auto-detection]",
                metavar="integer")
)

opt_parser = OptionParser(option_list = options)
args = parse_args(opt_parser)

# check validity of options ==========
if (is.null(args$d) | is.null(args$w) | is.null(args$o)){
    print_help(opt_parser)
    stop("All three arguments must be supplied.", call.=FALSE)
}

# Import and preprocess data ###############
d.raw <- read.delim(args$d, stringsAsFactors = FALSE)  # import all distances
samples <- unique(d.raw$sample)

# import assignment of weights of difference sources
dw <- read.delim(args$w, stringsAsFactors = FALSE)
dw <- dw[order(dw$weight, decreasing = FALSE), ]
weights <- dw$weight
names(weights) <- dw$source  # weights becomes a named and decreasing sorted integer vector
remove(dw)

# split the data frame by samples for parallel computing
d.sample <- lapply(samples, function(s) subset(d.raw, sample == s))
names(d.sample) <- samples  # lapply does not return a named list

# Filter repetitive distance measurements based on their weights for every sample ###############
distFilter <- function(d.high, d.all, lower.level) {
    # recursively filters low-weight measurements
    if (nrow(d.high) > 0) {
        if (lower.level > 0) {  # >= 1
            d.low <- d.all[[lower.level]]  # It may have zero row as the subset function does not return a NULL.

            # look for the first non-empty data frame of a lower level
            while (nrow(d.low) == 0 & lower.level > 1) {
                lower.level <- lower.level - 1
                d.low <- d.all[[lower.level]]  # keep searching until the data frame at a level is not empty
            }

            if (lower.level == 1 & nrow(d.low) == 0) {  # no other non-empty data frames exist below the level of d.high
                df <- d.high
            } else {  # lower.level >= 1 & nrow(d.low) > 0, etc.
                # do the actual merging stuff: removing query pairs in d.low when they are already present in d.high
                low.excl <- rep(FALSE, times = nrow(d.low))
                for (i in 1 : nrow(d.high)) {
                    a1 <- d.high$query1[i]
                    a2 <- d.high$query2[i]
                    # It is not necessary to try the other way because the query sequences are the same for every source in the same sample by the analysis design.
                    low.excl <- (d.low$query1 == a1 & d.low$query2 == a2) | low.excl
                }
                low.incl <- !low.excl

                if (sum(low.incl) > 0) {  # if d.low has some rows to contribute to d.high
                    d.low.filtered <- d.low[low.incl, ]
                    df <- distFilter(rbind(d.high, d.low.filtered), d.all, lower.level - 1)
                } else {  # All rows of d.low are found in d.high
                    df <- distFilter(d.high, d.all, lower.level - 1)
                }
            }

        } else {  # lower.level = 0; reaches the bottom of recursions
            df <- d.high
        }
    } else {  # nrow(d.high) = 0
        while (nrow(d.high) == 0 & lower.level > 0) {
            d.high <- d.all[[lower.level]]
            lower.level <- lower.level - 1
        }

        if (nrow(d.high) == 0 & lower.level == 0) {
            print("Error: empty distance measurements found for a sample, please check your input distance table.")
            stop()
        } else {
            df <- distFilter(d.high, d.all, lower.level)  # low.level was deducted by one in the while loop.
        }
    }

    return(df)
}

filterDistsByWeights <- function(ds, sources, source.num) {
    # A wrapper function that iteratively runs distFilter
    # [parameters] ds: an element of the list d.sample; w: weights
    d.sources <- lapply(sources, function(s) subset(ds, source == s))  # split by sources; some elements may be an empty data frame (nrow = 0, ncol > 0)
    names(d.sources) <- sources
    df <- distFilter(d.sources[[source.num]], d.sources, source.num - 1)  # source.num corresponds to the source of the highest weight

    return(df)
}

# determine how many cores will be used
if (args$c > 0) {
    np <- args$c
} else {
    np <- detectCores(logical = TRUE)  # detect the number of processors
}

cat(np, "cores are used for filtering distance measurements.\n", sep = " ")

# run the cluster
clr <- makeCluster(np)  # Initiate the cluster
clusterExport(clr, "distFilter", envir = environment())  # publishes variables to the cluster
d.filtered <- parLapply(clr, d.sample, filterDistsByWeights, names(weights), length(weights))  # process distances for every sample
stopCluster(clr)

# bind rows into a data frame
d <- do.call(rbind, d.filtered)
rownames(d) <- NULL

# save the result
write.table(d, file = args$o, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

print("Done!")
