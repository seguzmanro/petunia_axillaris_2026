#!/usr/bin/env Rscript

# Function to load or install required packages
ver_load_packages <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE, quietly = TRUE))
  need <- libs[req == FALSE]
  if(length(need) > 0) {
    install.packages(need, repos = "http://cran.us.r-project.org")
    lapply(need, require, character.only = TRUE, quietly = TRUE)
  }
}

ver_load_packages(c("argparse", "conStruct", "fields"))

parser <- ArgumentParser(description="Run a standalone conStruct model")
parser$add_argument('--freqs', required=TRUE, help='Path to frequencies CSV file')
parser$add_argument('--coords', required=TRUE, help='Path to coordinates CSV file')
parser$add_argument('--prefix', required=TRUE, help='Output prefix')
parser$add_argument('--K', type="integer", required=TRUE, help='K value to run')
parser$add_argument('--model_type', type="character", required=TRUE, choices=c('sp', 'nsp'), help='Model type (sp or nsp)')
parser$add_argument('--n_iter', type="integer", default=10000000, help='Number of MCMC iterations')

args <- parser$parse_args()

# Load Allele Frequencies
allele.frequencies <- as.matrix(read.csv(args$freqs, row.names = 1, check.names = FALSE))

# Load Coordinates
coords_raw <- read.csv(args$coords, row.names = 1, header = FALSE)
if (is.character(coords_raw[1, 1])) {
    coords_raw <- read.csv(args$coords, row.names = 1, header = TRUE)
    coords <- as.matrix(coords_raw[, 1:2])
} else {
    coords <- as.matrix(coords_raw)[, c(2, 1)]
}

# Calculate Geographic Distances
geoDist <- fields::rdist.earth(x1 = coords)
row.names(geoDist) <- row.names(coords)
colnames(geoDist) <- row.names(coords)

# Check if population order matches
if (!all(row.names(allele.frequencies) == row.names(coords))) {
    stop("ERROR: Population names/order in frequencies and coordinates do not match!")
}

cat(sprintf("Running %s conStruct model for K=%d...\n", args$model_type, args$K))

is_spatial <- (args$model_type == 'sp')

construct_run <- conStruct::conStruct(
    spatial = is_spatial,
    K = args$K,
    freqs = allele.frequencies,
    geoDist = geoDist,
    coords = coords,
    prefix = args$prefix,
    n.iter = args$n_iter,
    make.figs = TRUE,
    save.files = TRUE
)

cat("Run complete.\n")
