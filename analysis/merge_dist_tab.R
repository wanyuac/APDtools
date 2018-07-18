# This script merges a table of shortest-path distances (D) and a table of
# reference distances (Dr) into a single one. It also calculates the difference
# (D - Dr) and determines whether D is accurate under a given vector of thresholds.
#
# Parameters:
#   queries: filters distances for query sequences before sampling
#
# Example command line:
#   Rscript merge_dist_tab --dr Dr.tsv --d D.tsv --output sample1__Dm.tsv --cores 8
#
# Expect Dr to be produced by calc_dr.R in the ref_dist directory and expect D to be the output of
# compile_dists.py or parse_Bandage_output.sh.
#
# The sampling rates are ignored when the loci option is set.
#
# Copyright 2017 - 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 17-18 Dec 2017; the latest edition: 24 Mar 2018

# Read arguments ###############
library(optparse)
library(data.table)
library(parallel)

options = list(
    # input files
    make_option("--dr", type = "character", action = "store", default = NULL,
                help = "Path to the table of reference distances", metavar="character"),
    make_option("--d", type = "character", action = "store", default = NULL,
                help = "Path to the table of distance measurements", metavar = "character"),

    # overall filter
    make_option("--qs", type = "character", action = "store", default = "",
                help = "Comma-delimted names of query sequences whose distances will be merged."),

    # filters for shortest-path distances
    make_option("--max_d", type = "numeric", action = "store", default = 300e3,
                help = "Maximal shortest-path distance (bp) to be filtered for (before sampling). Set a negative value to turn the filter off. [default = %default]"),
    make_option("--max_nodes", type = "integer", action = "store", default = 10,
                help = "Maximal node number in distance paths to be filtered for before sampling. [default = %default] (<=0: turn off)"),
    make_option("--sampling_rate_d", type = "numeric", action = "store", default = 1,
                help = "Proportion (between 0 and 1) or number (>1) to sample from rows of the distance measurements. [default = %default] (turn off)"),

    # filters for reference distances
    make_option("--sampling_rate_dr", type = "numeric", action = "store", default = 1,
                help = "Proportion (between 0 and 1) or number (>1) to sample from rows of reference distances. [default = %default] (turn off)"),

    # executive arguments
    make_option("--output", type = "character", action = "store", default = "Dm.tsv",
                help = "Output file name and its path [default = %default]", metavar = "character"),
    make_option("--cores", type = "integer", action = "store", default = 4,
                help = "Number of processors to be used [default = %default]", metavar = "integer"),
    make_option("--saveRDS", type = "logical", action = "store_true", default = FALSE,
                help = "Turn it on to save the merged table as an RDS file instead of a tab-delimited file",
                metavar = "logical")
)

opt.parser = OptionParser(option_list = options)
args = parse_args(opt.parser)

# Function definitions ###############
keepQueries <- function(df, qs) {
    # Select rows of the data frame df by query1 and query2.
    # To increase the speed, this function takes two steps for filtering rather
    # than use the logic (query1 %in% qs) & (query2 %in% qs) because df usually
    # has too many rows.
    df <- subset(df, query1 %in% qs)  # reduce the number of rows for comparison in the next step
    df <- subset(df, query2 %in% qs)

    return(df)
}

sampleDistances <- function(df, rate = 1, suffix = "D", filename_prefix) {
    # Sample rows from df
    if (rate > 0) {
        # validity check
        n <- nrow(df)
        if (n > .Machine$integer.max) {  # a limit of the sample.int function
            print("Warning: the range of sampling does not cover all samples as the number of rows in df exceeds the machine's maximal integer.")
        }

        # draw indices for sampling
        if (rate < 1) {
            print(paste0("Sampling ", rate * 100, "% of rows from ", suffix, "."))
            indices <- sample.int(n = n, size = max(ceiling(n * rate), 1), replace = FALSE)  # At least there is one row gets selected.
        } else if (rate > 1) {
            print(paste0("Sampling ", rate, " rows from ", suffix, "."))
            indices <- sample.int(n = n, size = min(rate, n), replace = FALSE)  # Notice rate may equal n.
        } else {  # rate = 1
            print(paste("No sampling of", suffix,
                        "is performed as the rate equals one.", sep = " "))
        }

        # extract selected rows
        if (rate == 1) {
            df_sub <- df
        } else {
            df_sub <- df[sort(indices, decreasing = FALSE), ]
            print(paste("Saving the sampled data frame of", length(indices), "rows.", sep = " "))  # the number of extracted rows
            write.table(x = df_sub, file = paste(filename_prefix, suffix, "sampled.tsv", sep = "__"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        }
    } else {
        print("Warning: distances are not sampled because the sampling rate is non-positive.")
        df_sub <- df
    }
    rownames(df_sub) <- NULL

    return(df_sub)
}

getFilenamePrefix <- function(f) {
    # E.g. MG1655__Dm.tsv -> MG1655_Dm
    fields <- strsplit(x = f, split = ".", fixed = TRUE)[[1]]  # c("MG1655_Dm", "tsv")
    len <- length(fields)  # e.g., c("MG1655", "GenBank", "tsv"), where the last item is the filename extension.
    if (len > 1) {
        prefix <- paste(fields[-len], collapse = ".")
    } else {  # There is no "." in f.
        prefix <- fields
    }

    return(prefix)
}

mergeD2Dr <- function(r) {
    # Merge D (in the global environment) into Dr row-by-row. The argument r is a named vector from Dr.
    # D is large, hence do not pass D as an argument into each instance of this function.
    q1 <- r[["query1"]]
    q2 <- r[["query2"]]

    # First, look for combinations that are exactly the same.
    where <- which((D$query1 == q1) & (D$query2 == q2))
    if (length(where) > 0) {  # when there is a hit (cannot be more due to the algorithm for distance measurements)
        row_D <- D[where, ]  # get content of that row
        r_new <- data.frame(query1 = q1, query2 = q2, d = as.numeric(row_D$d),
                            node_number = row_D$node_number, dr = as.numeric(r[["dr"]]),
                            stringsAsFactors = FALSE)
    } else {  # search for the inverse order
        where <- which((D$query1 == q2) & (D$query2 == q1))
        if (length(where) > 0) {
            row_D <- D[where, ]
            r_new <- data.frame(query1 = q1, query2 = q2, d = as.numeric(row_D$d),
                                node_number = row_D$node_number, dr = as.numeric(r[["dr"]]),
                                stringsAsFactors = FALSE)  # keep the order defined in Dr
        } else {  # unique query combinations in Dr
            #r_new <- data.frame(query1 = q1, query2 = q2, d = NA,
            #                    node_number = NA, dr = as.numeric(r[["dr"]]),
			#					stringsAsFactors = FALSE)
            r_new <- NULL  # Do not add this row to the final data frame if we are not interested in missing measurability in de novo assemblies.
        }
    }

    return(r_new)  # empty or with a single row
}

isEmpty <- function(df, id = "D") {
    # Check if a data frame has no row.
    if (nrow(df) == 0) {
        stop(paste("Error:", id, "lost all rows after the filtering.", sep = " "))
    }
}

# Main program ###############
root <- getFilenamePrefix(args$output)

flt_qs <- args$qs != ""
if (flt_qs) {
    qs <- strsplit(x = args$qs, split = ",", fixed = TRUE)[[1]]  # a character vector
}

# Read inputs ===============
# Read reference distances from the output tables of calc_dr.R
print(paste("Reading", args$dr, sep = " "))
Dr <- fread(input = args$dr, sep = "\t", header = TRUE, data.table = FALSE,
            stringsAsFactors = FALSE, col.names = c("contig", "query1", "query2", "dr"))  # calc_dr.R prints four columns in the result.
Dr <- as.data.frame(Dr[, c("query1", "query2", "dr")], stringsAsFactors = FALSE)  # do not need contig information; convert data.table::data.frame to a standard data frame

# Filter and sampling of reference distances ===============
if (flt_qs) {
    print("Filtering reference distances for specific query sequences.")
    Dr <- keepQueries(Dr, qs)
    isEmpty(Dr, "Dr")
}
if (args$sampling_rate_dr != 1) {
    Dr <- sampleDistances(df = Dr, rate = args$sampling_rate_dr, suffix = "Dr",
                          filename_prefix = root)  # sample rows from Dr when it is specified to reduce the amount of time
} else {
    print("No sampling is performed as the rate equals one.")
}

# Read the output of compile_dists.py: it may be time-consuming as the input file is usually huge
print(paste("Reading", args$d, sep = " "))
D <- fread(input = args$d, sep = "\t", header = TRUE, data.table = FALSE,
           stringsAsFactors = FALSE)
D <- D[, c("query1", "query2", "distance", "node_number")]  # drop six columns that are not used so far (including the sample name column)
names(D)[3] <- "d"  # distance -> d
D <- as.data.frame(D, stringsAsFactors = FALSE)

# Apply filters for shortest-path distances
if (args$max_d >= 0) {
    print(paste0("Filtering D for a maximal distance of ", args$max_d, "bp."))
    D <- subset(D, d <= args$max_d)
    isEmpty(D, "D")
    print(paste("There are", nrow(D), "rows passed the distance filter.", sep = " "))
}
if (args$max_nodes > 0) {
    print(paste0("Filtering D for a maximum node number of ", args$max_nodes, "."))
    D <- subset(D, node_number <= args$max_nodes)
    isEmpty(D, "D")
    print(paste("There are", nrow(D), "rows passed the node-number filter.", sep = " "))
}
if (flt_qs) {  # Put this step to the last to save time as it involves the highest time complexity.
    print("Filtering shortest-path distances for specific query sequences.")
    D <- keepQueries(D, qs)
    isEmpty(D, "D")
}

# Sampling rows from the remaining data frame
if (args$sampling_rate_d != 1) {  # This is faster than calling sampleDistances when rate = 1.
    D <- sampleDistances(df = D, rate = args$sampling_rate_d, suffix = "D", filename_prefix = root)
} else {
    print("No sampling is performed as the rate equals one.")
}

# Merge D into Dr in accordance with query1 and query2 combinations in Dr ===============
# D may contain unique combinations of query1 and query2 because Dr does not have distances measured for
# queries that are not co-localised in the same molecule; Dr may also have unique query1-query2 combinations
# because there may be some distances unmeasurable in the assembly graph. Accordingly, this step ignores
# unique combinations of queries in D assuming Dr reveals all possible all-to-all distances (i.e. not include
# distances between loci on different DNA molecules). Some users may need to consider those unique combinations
# of queryies in D as well.
# In summary, the accuracy evaluated in this script is the conditional probability that
# a shortest-path distance is accurate when the corresponding pair of alleles are
# co-localised in a DNA molecule.
print(paste0(args$output, ": merging D into Dr using ", args$cores, " cores."))
cl <- makeCluster(min(args$cores, detectCores(logical = TRUE)))
clusterExport(cl = cl, varlist = "D", envir = environment())  # make variables accessible to different cores
Dm <- parApply(cl = cl, X = Dr, MARGIN = 1, mergeD2Dr)  # Merge D into Dr, ignoring query pairs that are not present in Dr.
stopCluster(cl)
Dm <- as.data.frame(rbindlist(Dm), stringsAsFactors = FALSE)
n <- nrow(Dm)
if (n > 0) {
    print(paste0(args$output, ": the merge step is completed and ", n, " rows are obtained."))
    Dm$error <- Dm$d - Dm$dr  # compute absolute errors; Any NA remains NA in the error column.

    # Store result ###############
    print(paste0(args$output, ": saving the result"))
    if (args$saveRDS) {
        saveRDS(Dm, file = args$output)
    } else {
        write.table(Dm, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)
    }
} else {
    print("The merge step is completed but there is no row mergable.")
}

print(paste(args$output, "has been produced successfully.", sep = " "))
