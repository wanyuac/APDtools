# Shared functions between accuracy_vs_dist.R and accuracy_vs_nodes.R
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 5 Jan 2018; the latest edition: 12 Jan 2018

getStrainNames <- function(filenames) {
    # Extract a strain name from every filename assuming the aforementioned format of filenames.
    # It returns a named vector with strain names as keys and file paths as values.
    filenames <- filenames[filenames != ""]  # For safty, even though the strsplit function does not return an empty element at the end of its return vector.
    strains <- vector(mode = "character", length = 0)
    for (f in filenames) {
        base_name <- strsplit(x = f, split = "/", fixed = TRUE)[[1]]
        base_name <- base_name[length(base_name)]  # get the last element
        strain_name <- strsplit(x = base_name, split = "__", fixed = TRUE)[[1]][1]
        strains[strain_name] <- f
    }

    return(strains)
}

getErrorTolerance <- function(errors) {
    max_errors <- as.numeric(strsplit(x = args$error_tol, split = ",", fixed = TRUE)[[1]]) * 1000  # error tolerance (bp)
    min_error_tol <- which(max_errors < 0)
    if (length(min_error_tol) > 0) {
        print("Warning: deleting negative cut-offs for error tolerance.")
        max_errors <- max_errors[-min_error_tol]
        if (length(max_errors) == 0) {
            stop("Error: cut-offs for errors must include non-negative values.")
        }
    }
    if (min(max_errors) > 0) {
        print("Warning: adding zero into the vector of tolerance cut-offs for exact measurements.")
        max_errors <- c(0, max_errors)  # Then the vector of error tolerance starts from zero.
    }

    return(max_errors)
}
