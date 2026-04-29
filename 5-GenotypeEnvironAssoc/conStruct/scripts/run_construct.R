#!/usr/bin/env Rscript

load_required_packages <- function(packages, lib_loc = NULL) {
  missing <- packages[!vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE,
    lib.loc = lib_loc
  )]

  if (length(missing) > 0) {
    stop(
      sprintf(
        "Missing required R packages: %s%s",
        paste(missing, collapse = ", "),
        if (!is.null(lib_loc) && nzchar(lib_loc)) {
          sprintf(" (searched in %s)", lib_loc)
        } else {
          ""
        }
      ),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    for (pkg in packages) {
      library(pkg, character.only = TRUE, lib.loc = lib_loc)
    }
  })
}


extract_cli_value <- function(flag_name) {
  args <- commandArgs(trailingOnly = TRUE)
  flag_idx <- match(flag_name, args)

  if (is.na(flag_idx)) {
    return(NULL)
  }

  value_idx <- flag_idx + 1
  if (value_idx > length(args)) {
    return(NULL)
  }

  args[[value_idx]]
}


cli_r_lib <- extract_cli_value("--r_lib")
load_required_packages(c("argparse", "conStruct", "fields"), lib_loc = cli_r_lib)

parser <- argparse::ArgumentParser(
  description = "Run a single conStruct analysis from an allele-frequency matrix and population coordinates."
)

parser$add_argument("--freqs", type = "character", required = TRUE, help = "Allele-frequency CSV")
parser$add_argument("--coords", type = "character", required = TRUE, help = "Population coordinate CSV")
parser$add_argument("--out_prefix", type = "character", required = TRUE, help = "Prefix for conStruct outputs")
parser$add_argument("--fit_out", type = "character", required = TRUE, help = "Output .rds file")
parser$add_argument("--done_file", type = "character", required = TRUE, help = "Completion marker file")
parser$add_argument("--k", type = "integer", required = TRUE, help = "Number of layers (K)")
parser$add_argument("--n_iter", type = "integer", default = 1000000, help = "Number of MCMC iterations")
parser$add_argument("--n_chains", type = "integer", default = 1, help = "Number of MCMC chains")
parser$add_argument("--coords_pop_col", type = "character", default = "Pop", help = "Population column in the coordinate table")
parser$add_argument("--coords_lon_col", type = "character", default = "lon", help = "Longitude column in the coordinate table")
parser$add_argument("--coords_lat_col", type = "character", default = "lat", help = "Latitude column in the coordinate table")
parser$add_argument("--spatial", action = "store_true", default = FALSE, help = "Run the spatial conStruct model")
parser$add_argument("--make_figs", action = "store_true", default = FALSE, help = "Allow conStruct to create built-in plots")
parser$add_argument("--save_files", action = "store_true", default = FALSE, help = "Allow conStruct to write its native output files")
parser$add_argument("--r_lib", type = "character", default = NULL, help = "Optional custom R library path")

args <- parser$parse_args()


read_frequency_matrix <- function(freq_path) {
  freq_df <- read.csv(freq_path, row.names = 1, check.names = FALSE)
  freq_mat <- as.matrix(freq_df)
  storage.mode(freq_mat) <- "double"
  freq_mat
}


prepare_coordinates <- function(coords_path, pop_col, lon_col, lat_col, pop_order) {
  coords_df <- read.csv(coords_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_cols <- c(pop_col, lon_col, lat_col)
  missing_cols <- setdiff(required_cols, colnames(coords_df))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Coordinate file %s is missing required columns: %s",
        coords_path,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coords_df <- coords_df[, required_cols]
  colnames(coords_df) <- c("population", "longitude", "latitude")
  coords_df$population <- as.character(coords_df$population)
  rownames(coords_df) <- coords_df$population

  missing_pops <- setdiff(pop_order, coords_df$population)
  if (length(missing_pops) > 0) {
    stop(
      sprintf(
        "Coordinate file %s is missing populations present in the allele-frequency matrix: %s",
        coords_path,
        paste(missing_pops, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coords_df <- coords_df[pop_order, c("longitude", "latitude"), drop = FALSE]
  coords_mat <- as.matrix(coords_df)
  storage.mode(coords_mat) <- "double"
  rownames(coords_mat) <- pop_order
  coords_mat
}


fit_dir <- dirname(args$fit_out)
if (!dir.exists(fit_dir)) {
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
}

out_dir <- dirname(args$out_prefix)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

freqs_mat <- read_frequency_matrix(args$freqs)
coords_mat <- prepare_coordinates(
  coords_path = args$coords,
  pop_col = args$coords_pop_col,
  lon_col = args$coords_lon_col,
  lat_col = args$coords_lat_col,
  pop_order = rownames(freqs_mat)
)

coords_used_path <- paste0(args$out_prefix, "_coords_used.csv")
write.csv(
  data.frame(
    population = rownames(coords_mat),
    longitude = coords_mat[, 1],
    latitude = coords_mat[, 2],
    row.names = NULL
  ),
  file = coords_used_path,
  row.names = FALSE,
  quote = FALSE
)

geo_dist <- NULL
if (isTRUE(args$spatial)) {
  geo_dist <- fields::rdist.earth(x1 = coords_mat)
  rownames(geo_dist) <- rownames(coords_mat)
  colnames(geo_dist) <- rownames(coords_mat)
}

construct_fit <- conStruct::conStruct(
  spatial = isTRUE(args$spatial),
  K = args$k,
  freqs = freqs_mat,
  geoDist = geo_dist,
  coords = coords_mat,
  prefix = args$out_prefix,
  n.chains = args$n_chains,
  n.iter = args$n_iter,
  make.figs = isTRUE(args$make_figs),
  save.files = isTRUE(args$save_files)
)

saveRDS(construct_fit, file = args$fit_out)

writeLines(
  c(
    sprintf("status\tcomplete"),
    sprintf("k\t%s", args$k),
    sprintf("spatial\t%s", isTRUE(args$spatial)),
    sprintf("freqs\t%s", normalizePath(args$freqs, winslash = "/", mustWork = TRUE)),
    sprintf("coords\t%s", normalizePath(args$coords, winslash = "/", mustWork = TRUE)),
    sprintf("fit_out\t%s", normalizePath(args$fit_out, winslash = "/", mustWork = FALSE)),
    sprintf("out_prefix\t%s", normalizePath(args$out_prefix, winslash = "/", mustWork = FALSE))
  ),
  con = args$done_file
)
