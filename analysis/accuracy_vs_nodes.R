# Draw a line graph to show accuracy of distance measurements (within 300 kb) as a function of node numbers (1..5).
#
# Input files are tab-delimited and must follow the same output format defined by the script merge_dist_tab.R.
# Expect their filenames to follow the format: [strain name]__suffix.tsv or ../../[strain name]__suffix.tsv.
#
# Example command line:
#   Rscript accuracy_vs_nodes.R --ds $(ls -1 *__Dm.tsv | tr '\n' ',')
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 5 Jan 2018; the latest version: 6 Jan 2018

# Load dependencies ###############
library(optparse)
library(ggplot2)
library(grid)
library(gridExtra)
library(parallel)
library(rlist)  # to use the function list.append
library(reshape2)

# Constants ###############
N_MAX <- 10  # only retain distances measured in paths of no more than five nodes
STRAIN_NMAX <- 12  # Currently, this script only handles a maximum of 12 samples.

# Obtain arguments ###############
options = list(
    # accuracy arguments
    make_option("--ds", dest = "ds", type = "character", default = NULL,
                help = "Output file names of merge_dist_tab.R, comma-delimited", metavar = "CHARACTER"),
    make_option("--d_breaks", dest = "d_breaks", type = "character", default = "1,5,10,50,100,200,250,300",
                help = "Breakpoints (kb) of distance measurements, comma-delimited [default = %default]", metavar = "CHARACTER"),
    make_option("--n_max", dest = "n_max", type = "integer", default = N_MAX,
                help = "Maximal node numbers to be accepted [default = %default]", metavar = "INTEGER"),
    make_option("--error_tol", dest = "error_tol", type = "character", default = "0,0.5,1,1.5,2,2.5",
                help = "Comma-delimited vector of error tolerance (kb) [default = %default]",
                metavar = "CHARACTER"),

    # settings for the output image
    make_option("--img", dest = "img", type = "character", default = "accuracy_vs_nodeNum.png",
                help = "Name (and its path) for the output image [default = %default]", metavar = "CHARACTER"),
    make_option("--width_per_panel", dest = "width_per_panel", type = "integer", default = 100,
                help = "Panel width in millimetres [default = %default]", metavar = "INTEGER"),
    make_option("--height_per_panel", dest = "height_per_panel", type = "integer", default = 40,
                help = "Panel height in millimetres [default = %default]", metavar = "INTEGER"),
    make_option("--res", dest = "res", type = "integer", default = 300,
                help = "Image resolution (ppi) [default = %default]", metavar = "INTEGER"),

    # installation directory
    make_option("--script_dir", dest = "script_dir", type = "character", default = ".",
                help = "Directory where func.R is located [default = %default] (no forward slash)"),

    # other options
    make_option("--cores", dest = "cores", type = "integer", action = "store", default = 4,
                help = "Number of processors to be used [default = %default]", metavar = "integer")
)

opt_parser = OptionParser(option_list = options)
args = parse_args(opt_parser)

# Functions ###############
# shared functions
source(paste(args$script_dir, "func.R", sep = "/"))  # get the script's directory and find func.R under the same directory

# private functions
importDistanceBreaks <- function(break_def) {
    # Import and perform sanity check on distance breakpoints
    d_breaks <- as.numeric(strsplit(x = break_def, split = ",", fixed = TRUE)[[1]])
    d_breaks <- sort(unique(d_breaks), decreasing = FALSE)
    where_invalid <- which(d_breaks < 0)
    if (length(where_invalid) > 0) {
        print("Warning: replacing negative distance breakpoints to zero.")
        d_breaks[where_invalid] <- 0
        d_breaks <- unique(d_breaks)
    }
    d_breaks <- d_breaks * 1000  # kb -> bp

    return(d_breaks)
}

determineAccuracy <- function(s, strains, n_max, d_breaks, error_tols) {
    require(data.table)  # for large input

    # import the merged distance table of the strain s and compute accuracy of distance measurements
    input_f <- strains[[s]]  # get the path to the input distance table
    print(paste(paste0("Processing the strain", s, "for statistics.", sep = " "), input_f, sep = ":"))
    output_f <- paste(s, "accuracy_nodes.tsv", sep = "__")  # [strain name]__accuracy_nodes.tsv
    if (file.exists(output_f)) {
        # Skip the time-consuming step for generating sliding-window statistics when it has been done
        print(paste("Reading statistics from the existing file:", output_f, sep = " "))
        ac <- fread(input = output_f, sep = "\t", header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)  # sliding-window statistics
    } else {
        # otherwise, process the data from scratch
        ds <- fread(input = input_f, sep = "\t", header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)
        ds <- ds[, c("query1", "query2", "d", "error", "node_number")]  # only use these five columns
        d_max <- max(d_breaks)
        ds <- subset(ds, d <= d_max)  # to remove extraordinarily large distances to reduce the amount of memory usage

        # determine accuracy for the current strain
        ac <- data.frame(Error_tol = integer(0), Node_num_max = integer(0), D_max = integer(0),
                         Accuracy = numeric(0), stringsAsFactors = FALSE)
        for (max_error in error_tols) {  # for each error tolerance level
            for (n_upper in 1 : n_max) {  # for each range of node numbers; n_max has been ensured to be no less than one.
                dn <- subset(ds, node_number <= n_upper)  # distances filtered for an upper bound of node numbers
                for (d_upper in d_breaks) {  # for each range of SPDs
                    if (d_upper < d_max) {
                        dnd <- subset(dn, d <= d_upper)  # Notice the distances have also been filtered for D_MAX.
                    } else {
                        dnd <- dn  # So do not need to repeat the previous subset process.
                    }
                    if (nrow(dnd) > 0) {  # The data frame dnd may be empty when d_upper is small.
                        a <- round(sum(abs(dnd$error) <= max_error) / nrow(dnd) * 100, digits = 2)  # then calculate the accuracy
                        ac <- rbind.data.frame(ac,
                                               data.frame(Error_tol = max_error,
                                                          Node_num_max = n_upper,
                                                          D_max = d_upper,
                                                          Accuracy = a,
                                                          stringsAsFactors = FALSE),
                                               stringsAsFactors = FALSE)
                    }  # do nothing otherwise
                }
            }
        }

        # save raw data of the current strain s
        write.table(ac, file = output_f, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    ac <- cbind.data.frame(Strain = rep(s, times = nrow(ac)), ac, stringsAsFactors = FALSE)  # Finally, append the strain name to the data frame ac.

    return(ac)
}

# Main program ###############
# extract strain names from file names
# Notice args$ds usually terminates by a comma due to the "tr" command, but luckily, the strsplit() function
# does not return an empty element at the end of its return vector.
strains <- getStrainNames(strsplit(x = args$ds, split = ",", fixed = TRUE)[[1]])  # returns a named vector of file paths
strain_num <- length(strains)
if (strain_num > STRAIN_NMAX) {
    print(paste("Warning: only the first", STRAIN_NMAX, "strains will be processed due to the limit of sample size.", sep = " "))
    strains <- strains[1 : STRAIN_NMAX]
    strain_num <- STRAIN_NMAX
}
strain_names <- names(strains)
print(paste0("Strains: ", paste(strain_names, collapse = ", "), "."))

# parse and check breakpoints of distances
d_breaks <- importDistanceBreaks(args$d_breaks)  # unit: bp
print(paste0("Breakpoints of the distances: ", paste(d_breaks, collapse = ", "), "."))

# check validity of node number
if (args$n_max < 1) {
    print(paste0("Warning: n_max must not be smaller than one. Resetting it to ", N_MAX, "."))
    args$n_max <- N_MAX
}

# determine valid error levels (bp)
error_cutoffs <- getErrorTolerance(args$error_tol)  # returns a vector of unnamed integers

# Calculate accuracy rate for each range of node numbers ===============
n_core <- min(args$cores, strain_num, detectCores())
cl <- makeCluster(n_core)
print(paste("Processing distance data of", strain_num, "strains using", n_core, "cores.", sep = " "))
clusterExport(cl = cl,
              varlist = list("strains", "error_cutoffs", "args", "d_breaks"),
              envir = environment())  # make variables accessible to different cores. "names" should be exported as well if you want to use names(strains).
ac_list <- parLapply(cl, strain_names, determineAccuracy,
                     strains, args$n_max, d_breaks, error_cutoffs)  # returns a list of data frames
stopCluster(cl)

# Reorganise ac_list into a named list ===============
ac <- vector(mode = "list", length = length(ac_list))
names(ac) <- strain_names
for (df in ac_list) {  # Columns: Strain, Error_tol, Node_num_max, D_max and Accuracy.
    s <- df$Strain[1]  # get the strain name
    df$Node_num_max <- factor(x = df$Node_num_max, levels = 1 : args$n_max,
                              labels = as.character(1 : args$n_max))
    ac[[s]] <- df[, -1]  # nrow(df) may equal zero when breakpoints in d_breaks are all small. Drop the Strain column.
}
rm(ac_list)  # release some memory

# Assign colours to different ranges of node numbers ===============
group_colours <- rainbow(n = args$n_max)
print("Saving the colour assignment for node counts.")
write.table(x = data.frame(N_max = 1 : args$n_max, Colour_RGBA = group_colours, stringsAsFactors = FALSE),
            file = paste0(args$img, "__colours.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

# Create ggplot objects ===============
panels <- list()  # ggplot objects
max_x <- max(d_breaks)

for (s in names(ac)) {  # go through strain names
    print(paste("Generating plot panels for the strain", s, sep = " "))
    for (error_level in error_cutoffs) {
        df <- subset(ac[[s]], Error_tol == error_level)  # accuracy under the current error level
        df <- df[, c("D_max", "Accuracy", "Node_num_max")]
        p <- ggplot(data = df) +
            geom_line(mapping = aes(x = D_max, y = Accuracy, group = Node_num_max, colour = Node_num_max)) +
            labs(title = paste0(s, ": +/-", round(error_level / 1000, digits = 1), " kb"),
                 x = "Max. of SPD (kb)", y = "Accuracy (%)") +
            scale_x_continuous(limits = c(0, max_x), breaks = c(0, d_breaks),
                               labels = as.character(c(0, round(d_breaks / 1000, digits = 1)))) +
            scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10),
                               labels = c("0", "", "20", "", "40", "", "60", "", "80", "", "100")) +
            scale_color_manual(breaks = as.character(1 : args$n_max), values = group_colours) + theme_bw() +
            theme(legend.position = "none", plot.title = element_text(size = 8, face = "bold"),
                  axis.text.x = element_text(size = 6, angle = 45, colour = "black"),
                  axis.text.y = element_text(size = 6, colour = "black"),
                  axis.title = element_text(size = 8, face = "bold"))
        panels <- rlist::list.append(panels, p)
    }

}

# Finally, draw ggplot objects into the output figure ===============
panel_columns <- length(error_cutoffs)
png(filename = args$img, units = "mm", res = args$res,
    width = args$width_per_panel * panel_columns,
    height = args$height_per_panel * strain_num)
par(oma = rep(0.1, times = 4), mar = rep(0.1, times = 4))
do.call("grid.arrange", c(panels, ncol = panel_columns))
dev.off()
