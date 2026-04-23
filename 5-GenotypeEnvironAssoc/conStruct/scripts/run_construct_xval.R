#!/usr/bin/env Rscript

# Function to load or install required packages
ver_load_packages <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE, quietly = TRUE))
  need <- libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need, repos = "http://cran.us.r-project.org")
    lapply(need, require, character.only = TRUE, quietly = TRUE)
  }
}

ver_load_packages(c("argparse", "conStruct", "fields", "doParallel", "foreach"))

parser <- ArgumentParser(description = "Run conStruct cross-validation")
parser$add_argument("--freqs", required = TRUE, help = "Path to frequencies CSV file")
parser$add_argument("--coords", required = TRUE, help = "Path to coordinates CSV file")
parser$add_argument("--prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--min_k", type = "integer", default = 1, help = "Minimum K to test")
parser$add_argument("--max_k", type = "integer", default = 8, help = "Maximum K to test")
parser$add_argument("--n_reps", type = "integer", default = 8, help = "Number of cross-validation replicates")
parser$add_argument("--train_prop", type = "double", default = 0.9, help = "Proportion of data to use for training")
parser$add_argument("--n_iter", type = "integer", default = 1000, help = "Number of MCMC iterations per replicate")
parser$add_argument("--threads", type = "integer", default = 1, help = "Number of threads")

args <- parser$parse_args()

# Load Allele Frequencies
allele.frequencies <- as.matrix(read.csv(args$freqs, row.names = 1, check.names = FALSE, header = TRUE))

# Load Coordinates
# We check if the first row contains characters to guess if there is a header.
coords_raw <- read.csv(args$coords, row.names = 1, header = FALSE)
if (is.character(coords_raw[1, 1]) || is.factor(coords_raw[1, 1])) {
  # It likely has a header. Read again with header=TRUE
  coords_raw <- read.csv(args$coords, row.names = 1, header = TRUE)
  coords <- as.matrix(coords_raw[, 1:2])
} else {
  # If no header, original script used columns 2 and 1 in that order.
  coords <- as.matrix(coords_raw)[, c(2, 1)]
}

# Ensure all populations in frequencies are present in coordinates
if (!all(row.names(allele.frequencies) %in% row.names(coords))) {
  missing_pops <- row.names(allele.frequencies)[!row.names(allele.frequencies) %in% row.names(coords)]
  stop(paste("ERROR: Some populations in frequencies are missing from coordinates:", paste(missing_pops, collapse=", ")))
}

# Re-order coordinates to exactly match the order of populations in allele frequencies
coords <- coords[row.names(allele.frequencies), ]

# Calculate Geographic Distances
geoDist <- fields::rdist.earth(x1 = coords)
row.names(geoDist) <- row.names(coords)
colnames(geoDist) <- row.names(coords)

# Final sanity check (should always pass now)
if (!all(row.names(allele.frequencies) == row.names(coords))) {
  stop("ERROR: Population names/order in frequencies and coordinates do not match!")
}

cat("Running conStruct cross-validation...\n")

# Run Cross Validation
xvals <- conStruct::x.validation(
  train.prop = args$train_prop,
  n.reps = args$n_reps,
  K = args$min_k:args$max_k,
  freqs = allele.frequencies,
  data.partitions = NULL,
  geoDist = geoDist,
  coords = coords,
  prefix = args$prefix,
  n.iter = args$n_iter,
  make.figs = FALSE,
  save.files = TRUE,
  parallel = TRUE,
  n.nodes = args$threads
)

cat("Cross-validation complete.\n")
