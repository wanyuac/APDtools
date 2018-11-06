# Draw a line graph to show accuracy of distance measurements (within 300 kb) as
# a function of the measurements
#
# For each pair of alleles:
# Accuracy = P(accurate distance measurements | a maximum of the measurements AND
# the alleles are co-localised in an actual genome).
#
# Input files are tab-delimited and must follow the same output format of the
# script merge_dist_tab.R. Expect their filenames to follow the format:
# [strain name]__suffix.tsv or ../../[strain name]__suffix.tsv.
#
# Example command line:
#   Rscript accuracy_vs_dist.R --ds $(ls -1 dists/*__Dm.tsv | tr '\n' ',')
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 20 Dec 2017, the latest edition: 6 Nov 2018

# Load dependencies ###############
library(optparse)
library(ggplot2)
library(grid)
library(gridExtra)
library(parallel)
library(reshape2)
library(data.table)

# Constants ###############
D_MAX <- 3e5  # only assess the measurements within 300 kb
IMG_NCOL <- 3  # number of columns in the output image
STRAIN_NMAX <- 12  # Currently, this script only handles a maximum of 12 samples.

# Obtain arguments ###############
options = list(
    # accuracy arguments
    make_option("--ds", dest = "ds", type = "character", default = NULL,
                help = "Comma-delimited output file names of merge_dist_tab.R",
                metavar = "CHARACTER"),
    make_option("--d_max", dest = "d_max", type = "integer", default = 50,
                help = "Maximal shrtest-path distances (kb) to be analysed [default = %default]",
                metavar = "INTEGER"),
    make_option("--error_tol", dest = "error_tol", type = "character", default = "0,0.5,1,1.5,2,2.5",
                help = "Comma-delimited vector of error tolerance in the unit of kb [default = %default]",
                metavar = "CHARACTER"),
    make_option("--window_size", dest = "window_size", type = "integer", default = 1000,
                help = "Size of sliding window in bp [default = %default]",
                metavar = "INTEGER"),

    # a filter for node numbers in every distance path
    make_option("--node_num", dest = "node_num", type = "integer", default = 1,
                help = "Maximal number of nodes per distance path allowed [default = %default]", metavar = "INTEGER"),
    make_option("--le", dest = "le", type = "logical", action = "store_true", default = FALSE,
                help = "Turn it on to select distances of node numbers less than or equal to the node_num", metavar = "logical"),

    # image settings
    make_option("--img", dest = "img", type = "character", default = "accuracy_vs_dist.png",
                help = "Name (and its path) for the output image [default = %default]", metavar = "CHARACTER"),
    make_option("--width", dest = "width", type = "integer", default = 150,
                help = "Image width in millimetres [default = %default]", metavar = "INTEGER"),
    make_option("--height", dest = "height", type = "integer", default = 220,
                help = "Image height in millimetres [default = %default]", metavar = "INTEGER"),
    make_option("--res", dest = "res", type = "integer", default = 300,
                help = "Image resolution (ppi) [default = %default]", metavar = "INTEGER"),

    # installation directory
    make_option("--script_dir", dest = "script_dir", type = "character", default = ".",
                help = "Directory where func.R is located [default = %default] (no forward slash)",
                metavar = "CHARACTER"),

    # other options
    make_option("--cores", dest = "cores", type = "integer", action = "store", default = 4,
                help = "Number of processors to be used [default = %default]", metavar = "integer")
)  # create a list of OptionParserOption instances to used as the option_list by the OptionParser function to create a new OptionParser instance

opt_parser = OptionParser(option_list = options)  # load the list of OptionParserOption instances to create an instance of the parser object
args = parse_args(opt_parser)  # parses command-line options using this OptionParser instance

# Functions ###############
# shared functions
source(paste(args$script_dir, "func.R", sep = "/"))  # get the script's directory and find func.R under the same directory

# private functions
determineAccuracy <- function(ds, max_errors) {
    # The accuracy evaluated in this script is the conditional probability that
    # a shortest-path distance is accurate when the corresponding pair of alleles
    # are co-localised in a DNA molecule.
    max_errors <- sort(max_errors, decreasing = FALSE)

    # count how many accuracy levels
    min_threshold <- max_errors[1]
    if (min_threshold < 0) {
        stop("Error: the maximum of errors must not be negative.")
    } else {
        if (min_threshold == 0) {
            ac_levels <- 0 : (length(max_errors) - 1)  # already include the "exact" level
        } else {
            ac_levels <- 0 : length(max_errors)  # add an additional column for exact measurements
            max_errors <- c(0, max_errors)
        }
    }

    # append columns for accuracy assertments
    # The accuracy is determined for the absolute value of errors.
    for (lv in ac_levels) {
        ds[[paste0("L", lv)]] <- !(abs(ds$error) > max_errors[lv + 1])  # append a new column to the data frame ds with names L0, L1, ...
    }

    return(list(tab = ds, levels = ac_levels))
}

calcStats <- function(ds, error_levels, w, d_max) {
    # Calculates summary statistics for sliding windows of a width w
    acc_col <- paste0("L", error_levels)  # column names for different accuracy levels
    d_max <- min(d_max, max(ds$d))  # actual maximum of distance measurements
    n_win <- d_max - w + 1  # number of windows
    win_ids <- 1 : n_win  # window ID, which equals the start coordinate of each window
    delta_win <- w - 1  # number of bases added to the current start base to get the upper bound of a window
    sw <- NULL  # initialise the result data frame

    # calculate accuracy in each window
    for (start in win_ids) {
        end <- start + delta_win
        midpoint <- (start + end) / 2  # the middle of each window, which will be used as coordinates on the X axis in the output plot
        dw <- subset(ds, d >= start & d <= end)  # distance measurements within the current window
        n <- nrow(dw)  # number of distance measurements within this window
        df <- data.frame(start = start, end = end, midpoint = midpoint, n = n, stringsAsFactors = FALSE)

        # append columns about accuracy to the temporary data frame df
        for (group in acc_col) {
            if (n > 0) {
                df[, group] <- round(sum(dw[, group]) / n, digits = 4) * 100
            } else {  # no measurements within the current window
                df[, group] <- NA
            }
        }
        sw <- rbind.data.frame(sw, df)
    }

    return(sw)
}

summariseSlidingWindows <- function(s, strains, max_errors, d_max, w, node_num, is_le) {
    # This is a function designed for parallel computing. w: window size
    # It takes as input data of a single strain.
    # Arguments:
    # w: window size; is_le: is less than or equal to the node_num

    input_f <- strains[[s]]
    print(paste(paste0("Processing the strain", s, "for statistics.", sep = " "),
                input_f, sep = ":"))
    output_f <- paste(s, "accuracy_dist_sw.tsv", sep = "__")  # [strain name]__accuracy_dist_sw.tsv
    if (file.exists(output_f)) { # Skip the time-consuming step for generating sliding-window statistics when it has been done
        print(paste("Reading statistics from the existing file:", output_f, sep = " "))
        sw <- fread(input = output_f, sep = "\t", header = TRUE, data.table = FALSE,
                    stringsAsFactors = FALSE)  # sliding-window statistics
    } else {  # process the data from scratch
        ds <- fread(input = input_f, sep = "\t", header = TRUE, data.table = FALSE,
                    stringsAsFactors = FALSE)
        ds <- ds[, c("query1", "query2", "d", "error", "node_number")]  # only use these five columns

        # selecting for node numbers
        if (is_le) {
            ds <- subset(ds, node_number <= node_num)
        } else {
            ds <- subset(ds, node_number == node_num)
        }

        # selecting distances for their sizes
        ds <- ds[, -5]  # only use the first four columns
        ds <- subset(ds, d <= d_max)  # drop shortest-path distances > d_max
        ds <- determineAccuracy(ds, max_errors)  # append columns to the original ds; This function returns a list.
        sw <- calcStats(ds = ds[["tab"]], error_levels = ds[["levels"]], w = w,
                        d_max = d_max)
        write.table(sw, file = output_f, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    sw <- cbind.data.frame(strain = rep(s, times = nrow(sw)), sw,
                           stringsAsFactors = FALSE)  # append the strain name for rbind.data.frame

    return(sw)
}

# Main program ###############
strains <- getStrainNames(strsplit(x = args$ds, split = ",", fixed = TRUE)[[1]])
strain_num <- length(strains)
if (strain_num > STRAIN_NMAX) {
    print(paste("Warning: only the first", STRAIN_NMAX, "strains will be processed due to the limit of sample size.", sep = " "))
    strains <- strains[1 : STRAIN_NMAX]
    strain_num <- STRAIN_NMAX
}

img_nrow <- ceiling(strain_num / IMG_NCOL)  # number of rows for panels in the image
d_max <- min(args$d_max * 1000, D_MAX)  # kb -> bp
print(paste("Shortest-path distances within", d_max, "bp will be analysed.", sep = " "))

if (args$node_num < 1) {
    print("Warning: the minimum of node_num must be at least one. Resetting this argument to one.")
    args$node_num <- 1
}

# determine valid error levels ===============
max_errors <- getErrorTolerance(args$error_tol)  # unit: bp

# Calculate statistics based on sliding windows ===============
n_core <- min(args$cores, strain_num, detectCores())
cl <- makeCluster(n_core)
print(paste("Processing distance data of", strain_num, "strains using", n_core, "cores.", sep = " "))
clusterExport(cl = cl,
              varlist = list("strains", "max_errors", "d_max", "args", "determineAccuracy", "calcStats"),
              envir = environment())  # make variables accessible to different cores
sw_cl <- parLapply(cl, names(strains), summariseSlidingWindows, strains, max_errors, d_max, args$window_size, args$node_num, args$le)
stopCluster(cl)

# Reorganise sw_cl into a named list ===============
sw <- list()  # sliding-window result for strains
for (tab in sw_cl) {
    s <- tab$strain[1]  # get the strain name
    sw[[s]] <- tab
}
rm(sw_cl)  # release some memory

# Create ggplot objects ===============
panels <- list()  # Do not use vector(mode = "list", length = length(sw)) as it will have three empty elements later.
for (s in names(sw)) {  # s: a strain name
    print(paste("Generating a plot for the strain", s, sep = " "))
    error_tol <- names(sw[[s]])
    error_tol <- error_tol[6 : length(error_tol)]  # chop off the first five columns, leaving only names of error levels L0, L1, ..., Ln
    stats <- reshape2::melt(data = sw[[s]], id.vars = "midpoint",
                            measure.vars = error_tol, variable.name = "error_tol",
                            value.name = "accuracy", na.rm = FALSE)
    group_colours <- rainbow(n = length(error_tol))
    x_breaks <- seq(0, d_max, length.out = 6)  # ticks dividing the X axis into five segments
    panels[[s]] <- ggplot(data = stats) +
        geom_line(mapping = aes(x = midpoint, y = accuracy, group = error_tol, colour = error_tol)) +
        labs(x = "Distance (kb)", y = "Accuracy (%)", title = s) +
        scale_x_continuous(limits = c(0, d_max), breaks = x_breaks,
                           labels = round(x_breaks / 1000, digits = 1)) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10),
                           labels = c("0", "", "20", "", "40", "", "60", "", "80", "", "100")) +
        scale_color_manual(breaks = error_tol, values = group_colours) + theme_bw() +
        theme(legend.position = "none", plot.title = element_text(size = 8, face = "bold"),
              axis.text.x = element_text(size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title = element_text(size = 8, face = "bold"))
}

# Draw ggplot objects into the output figure ===============
# Add the argument type = "cairo-png" to avoid the error: "unable to open connection to X11 display" on some systems.
# See https://stackoverflow.com/questions/24999983/r-unable-to-start-device-png-capabilities-has-true-for-png for details.
png(filename = args$img, width = args$width, height = args$height, units = "mm",
    res = args$res, type = "cairo-png")
par(oma = rep(0.1, times = 4), mar = rep(0.1, times = 4))
do.call("grid.arrange", c(panels, ncol = 3))
dev.off()

# Print colour assignment ===============
print("Printing colour assignment.")
write.table(x = data.frame(Max_error_bp = max_errors, Colour_RGBA = group_colours, stringsAsFactors = FALSE),
            file = paste0(args$img, "__colours.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
