#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(BEDASSLE)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(argparse)
  library(dplyr)
  library(coda)
  library(ggmcmc)
})

script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub('^--file=', '', script_arg[1])))
source(file.path(script_dir, 'bedassle_diagnostics.R'))

# Helper function to convert pairwise distance table to a distance matrix
dist_vec_to_matrix <- function(pairwise_df, var_name, pop_names) {
  mat <- matrix(0, nrow=length(pop_names), ncol=length(pop_names), dimnames=list(pop_names, pop_names))
  for (i in 1:nrow(pairwise_df)) {
    p1 <- as.character(pairwise_df$p1[i])
    p2 <- as.character(pairwise_df$p2[i])
    if (p1 %in% pop_names && p2 %in% pop_names) {
      val <- as.numeric(pairwise_df[[var_name]][i])
      mat[p1, p2] <- val
      mat[p2, p1] <- val
    }
  }
  return(mat)
}

# BEDASSLE requires strictly non-negative distances, with 0 on the diagonal.
# We shift the matrix so its minimum off-diagonal is 0, just in case the input was centered (scaled to have mean 0).
adjust_dist_matrix <- function(mat) {
  diag(mat) <- NA
  min_val <- min(mat, na.rm = TRUE)
  if (min_val < 0) {
    mat <- mat - min_val
  }
  diag(mat) <- 0
  return(mat)
}

parser <- ArgumentParser(description='Refactored BEDASSLE R script for Pop Genomics')

parser$add_argument('--vcf', type='character', help='VCF file (required for new MCMC runs)')
parser$add_argument('--popmap', type='character', help='Population map CSV (required for new MCMC runs)')
parser$add_argument('--env_dist', type='character', required=TRUE, help='CSV file with pairwise geographic and environmental distances')
parser$add_argument('--out_prefix', type='character', required=TRUE, help='Output prefix')

# MCMC Parameters
parser$add_argument('--ngen', type='integer', default=17000000, help='Total MCMC iterations, including burn-in')
parser$add_argument('--burnin', type='integer', default=4250000, help='Iterations to exclude before summaries and diagnostics')
parser$add_argument('--min_saved_states', type='integer', default=50000, help='Minimum post-burn-in saved states for new runs')
parser$add_argument('--ess_min', type='double', default=400, help='ESS screening threshold; flags do not prove convergence')
parser$add_argument('--seed', type='integer', default=20261006, help='Base seed; each predictor receives a distinct seed')
parser$add_argument('--summarize_only', action='store_true', help='Reprocess completed saved chains without running MCMC')
parser$add_argument('--summary_input_prefix', type='character', help='Saved-chain prefix for summary-only mode')
parser$add_argument('--no_plots', action='store_true', help='Write numeric diagnostics without PDF reports')
parser$add_argument('--printfreq', type='integer', default=10000, help='Print frequency')
parser$add_argument('--savefreq', type='integer', default=100000, help='Save frequency')
parser$add_argument('--samplefreq', type='integer', default=250, help='Sample frequency')
parser$add_argument('--delta', type='double', default=1e-20, help='Delta (small positive number)')
parser$add_argument('--aD_stp', type='double', default=0.075, help='aD step size')
parser$add_argument('--aE_stp', type='double', default=0.05, help='aE step size')
parser$add_argument('--a2_stp', type='double', default=0.025, help='a2 step size')
parser$add_argument('--phi_stp', type='double', default=30, help='phi step size')
parser$add_argument('--thetas_stp', type='double', default=0.2, help='thetas step size')
parser$add_argument('--mu_stp', type='double', default=0.35, help='mu step size')

parser$add_argument('--threads', type='integer', default=1, help='Number of threads for parallel MCMCs')

args <- parser$parse_args()

if (!is.finite(args$ess_min) || args$ess_min <= 0 || args$threads < 1 ||
    args$seed < 0 || as.double(args$seed) + 10000 > .Machine$integer.max) {
  stop('Invalid ESS threshold, thread count, or seed.')
}
proposal_settings <- unlist(args[c('delta', 'aD_stp', 'aE_stp', 'a2_stp',
                                   'phi_stp', 'thetas_stp', 'mu_stp')])
if (any(!is.finite(proposal_settings)) || any(proposal_settings <= 0) || args$printfreq < 1) {
  stop('Proposal scales, delta and print frequency must be positive.')
}

out_dir <- dirname(args$out_prefix)
if (out_dir != "." && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
out_dir <- normalizePath(out_dir)
args$out_prefix <- file.path(out_dir, basename(args$out_prefix))

# Model order and summary inputs are defined by the environmental-distance table.
env_dist <- read.csv(args$env_dist, check.names = FALSE)
if (!all(c('pop1', 'pop2', 'geog') %in% names(env_dist))) stop('Distance table lacks required columns.')
env_var_names <- setdiff(colnames(env_dist), c('pop1', 'pop2', 'p1', 'p2', 'geog'))
if (!length(env_var_names) || anyDuplicated(env_var_names)) stop('Invalid environmental predictors.')

diagnostics_dir <- paste0(args$out_prefix, '_mcmc_plots')
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), paste0(args$out_prefix, '_session_info.txt'))
saveRDS(args, paste0(args$out_prefix, '_postprocessing_settings.rds'))

write_bedassle_summaries <- function(input_prefix) {
  summaries <- lapply(env_var_names, function(variable) {
    path <- paste0(input_prefix, '_', variable, '_MCMC_output1.Robj')
    summarize_bedassle_file(path, variable, args$burnin,
      file.path(diagnostics_dir, variable), ess_min = args$ess_min,
      min_saved_states = args$min_saved_states, make_plots = !args$no_plots)
  })
  table <- do.call(rbind, lapply(summaries, `[[`, 'summary'))
  parameters <- do.call(rbind, lapply(summaries, `[[`, 'diagnostics'))
  # Write aggregate files only once every requested model has been processed.
  write.csv(parameters, paste0(args$out_prefix, '_parameter_diagnostics.csv'), row.names = FALSE)
  write.csv(table, paste0(args$out_prefix, '_BEDASSLE_RES_CI.csv'), row.names = FALSE)
  cat(sprintf('Saved posterior credible intervals and diagnostics: %s\n', args$out_prefix))
}

if (args$summarize_only) {
  input_prefix <- args$summary_input_prefix
  if (is.null(input_prefix)) stop('--summary_input_prefix is required in summary-only mode.')
  input_prefix <- file.path(normalizePath(dirname(input_prefix)), basename(input_prefix))
  if (identical(input_prefix, args$out_prefix)) stop('Use a different output prefix to preserve original summaries.')
  write_bedassle_summaries(input_prefix)
  quit(save = 'no', status = 0)
}

if (is.null(args$vcf) || is.null(args$popmap)) stop('--vcf and --popmap are required for MCMC runs.')
retained <- validate_bedassle_settings(args$ngen, args$burnin, args$samplefreq,
                                      args$savefreq, args$min_saved_states)
cat(sprintf('MCMC iterations: %d; burn-in: %d; retained states per model: %d\n',
            args$ngen, args$burnin, retained))
existing <- paste0(args$out_prefix, '_', env_var_names, '_MCMC_output1.Robj')
if (any(file.exists(existing))) stop('MCMC outputs already exist under this prefix. Use a new prefix or summary-only mode.')

run_metadata <- list(arguments = args,
  input_md5 = tools::md5sum(c(args$vcf, args$popmap, args$env_dist)),
  session_info = sessionInfo(), started_at = Sys.time(),
  diagnostic_scope = 'sampled covariance hyperparameters, beta, population phi, and aE/aD; locus-specific theta/mu trajectories are not saved by BEDASSLE')
saveRDS(run_metadata, paste0(args$out_prefix, '_run_metadata.rds'))

cat("Loading files...\n")
samples_info <- read.csv(args$popmap)
# First col: Indv, Second Col: Pop
colnames(samples_info)[1:2] <- c('indv', 'pop')
rownames(samples_info) <- samples_info$indv

loaded_vcf <- read.vcfR(args$vcf)
loaded_genind <- vcfR2genind(loaded_vcf)

# Retain only samples found in popmap
common_samples <- intersect(rownames(loaded_genind@tab), rownames(samples_info))
loaded_genind@tab <- loaded_genind@tab[common_samples, ]
samples_info <- samples_info[common_samples, ]

loaded_genind@pop <- as.factor(samples_info$pop)
loaded_genpop <- genind2genpop(loaded_genind, samples_info$pop)

# Convert to BEDASSLE format
bedassle_input <- as.matrix(loaded_genpop@tab)
del <- seq(2, ncol(bedassle_input), 2)
if (length(del) > 0) {
  bedassle_input <- bedassle_input[, -del] 
}

# Create matrix for sample sizes
samples_matrix_n <- matrix(nrow=nrow(bedassle_input), ncol=ncol(bedassle_input))
rownames(samples_matrix_n) <- rownames(bedassle_input)
size_sample <- table(loaded_genind@pop)
size_sample <- size_sample[rownames(samples_matrix_n)] # Retain specific order
size_sample_vector <- as.vector(size_sample)

for(i in 1:nrow(samples_matrix_n)){
  samples_matrix_n[i,] <- size_sample_vector[i]
}
samples_matrix_n <- samples_matrix_n * 2 # Adjust (diploid, loss of one allele)

cat("BEDASSLE input matrices prepared.\n")

# Parse Pairwise Distance Matrix (env_dist)
env_dist <- env_dist %>% rowwise() %>% mutate(
  p1 = min(as.character(pop1), as.character(pop2)),
  p2 = max(as.character(pop1), as.character(pop2))
) %>% ungroup()

pop_names_in_data <- rownames(bedassle_input)

# Prepare Geographic Distance Matrix D
geo_matrix <- dist_vec_to_matrix(env_dist, "geog", pop_names_in_data)
geo_matrix <- adjust_dist_matrix(geo_matrix)

# Critically, BEDASSLE's spatial covariance priors and user-defined step sizes 
# heavily depend on distances being statistically scaled. Since geographic distance is imported 
# natively in Km, we must divide the matrix by its own standard deviation to align it.
geo_matrix <- geo_matrix / sd(c(geo_matrix))

cat("Geographic distance matrix D constructed and standardized by SD.\n")

# Prepare Environmental Distance Matrices E
env_var_names <- setdiff(colnames(env_dist), c("pop1", "pop2", "p1", "p2", "geog"))
env_matrices <- list()

for (var in env_var_names) {
  e_mat <- dist_vec_to_matrix(env_dist, var, pop_names_in_data)
  e_mat <- adjust_dist_matrix(e_mat)
  env_matrices[[var]] <- e_mat
}
cat(sprintf("Prepared %d environmental variables for E.\n", length(env_matrices)))

# MCMC wrapper function
run_bedassle_mcmc <- function(E_matrix, var_name, index) {
  prefix_str <- paste0(args$out_prefix, "_", var_name, "_")
  old_dir <- getwd()
  log_connection <- file(paste0(prefix_str, 'chain.log'), open = 'wt')
  sink(log_connection)
  sink(log_connection, type = 'message')
  on.exit({
    sink(type = 'message')
    sink()
    close(log_connection)
    setwd(old_dir)
  }, add = TRUE)
  set.seed(args$seed + index)
  cat('Running BEDASSLE for variable:', var_name, '; seed:', args$seed + index, '\n')
  started <- Sys.time()
  
  # Note: BEDASSLE outputs to the working directory / directory argument
  res <- BEDASSLE::MCMC_BB(
    counts = bedassle_input,
    sample_sizes = samples_matrix_n,
    D = geo_matrix,
    E = E_matrix,
    k = nrow(bedassle_input), 
    loci = ncol(bedassle_input),
    delta = args$delta,
    aD_stp = args$aD_stp,
    aE_stp = args$aE_stp,
    a2_stp = args$a2_stp,
    phi_stp = args$phi_stp,
    thetas_stp = args$thetas_stp,
    mu_stp = args$mu_stp,
    ngen = args$ngen,
    printfreq = args$printfreq,
    savefreq = args$savefreq,
    samplefreq = args$samplefreq,
    directory = out_dir,
    prefix = basename(prefix_str),
    continue = FALSE,
    continuing.params = FALSE
  )
  cat('MCMC returned:', res, '\nElapsed seconds:',
      as.numeric(difftime(Sys.time(), started, units = 'secs')), '\n')
  load_bedassle_chain(paste0(prefix_str, 'MCMC_output1.Robj'), args$burnin)
  return(prefix_str)
}

num_threads <- min(args$threads, length(env_matrices))
cat(sprintf("Starting MCMC chains with %d parallel threads...\n", num_threads))
cl <- makeCluster(num_threads, type = 'FORK')
registerDoParallel(cl)

# BEDASSLE execution loop
prefixes <- tryCatch({
  foreach(iteration=seq_along(env_matrices), .errorhandling = 'stop') %dopar% {
    run_bedassle_mcmc(env_matrices[[iteration]], names(env_matrices)[iteration], iteration)
  }
}, finally = stopCluster(cl))
cat("All BEDASSLE MCMC runs completed.\n")

# Combine Results into a summary CSV
write_bedassle_summaries(args$out_prefix)
