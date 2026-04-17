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
})

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

parser$add_argument('--vcf', type='character', required=TRUE, help='VCF file')
parser$add_argument('--popmap', type='character', required=TRUE, help='Population map CSV')
parser$add_argument('--env_dist', type='character', required=TRUE, help='CSV file with pairwise geographic and environmental distances')
parser$add_argument('--out_prefix', type='character', required=TRUE, help='Output prefix')

# MCMC Parameters
parser$add_argument('--ngen', type='integer', default=2000000, help='Number of MCMC steps')
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

out_dir <- dirname(args$out_prefix)
if (out_dir != "." && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

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
env_dist <- read.csv(args$env_dist)
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
run_bedassle_mcmc <- function(E_matrix, var_name) {
  prefix_str <- paste0(args$out_prefix, "_", var_name, "_")
  cat("Running BEDASSLE for variable:", var_name, "\n")
  
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
  return(prefix_str)
}

num_threads <- min(args$threads, length(env_matrices))
cat(sprintf("Starting MCMC chains with %d parallel threads...\n", num_threads))
cl <- makeCluster(num_threads, type = 'FORK')
registerDoParallel(cl)

# BEDASSLE execution loop
prefixes <- foreach(iteration=1:length(env_matrices)) %dopar% {
  run_bedassle_mcmc(env_matrices[[iteration]], names(env_matrices)[iteration])
}

stopCluster(cl)
cat("All BEDASSLE MCMC runs completed.\n")

# Combine Results into a summary CSV
vars <- c()
var_means <- c()
var_CI <- c()

for (i in 1:length(env_matrices)) {
  var <- names(env_matrices)[i]
  # output structure: {directory}/{prefix}MCMC_output1.Robj
  robj_file <- file.path(out_dir, paste0(basename(args$out_prefix), "_", var, "_MCMC_output1.Robj"))
  
  if (file.exists(robj_file)) {
    load(robj_file) # Loads aE, aD vectors from MCMC footprint
    
    ratio <- aE / aD
    m_val <- round(mean(ratio, na.rm=TRUE), 4)
    if (length(ratio) > 1 && sd(ratio) > 0) {
      margin <- qt(0.975, df=(length(ratio)-1)) * sd(ratio) / sqrt(length(ratio))
      ci_str <- paste0('[', round(m_val - margin, 3), ' - ', round(m_val + margin, 3), ']')
    } else {
      # If fewer steps or 0 variance
      ci_str <- paste0('[', m_val, ' - ', m_val, ']')
    }
    
    vars <- c(vars, var)
    var_means <- c(var_means, m_val)
    var_CI <- c(var_CI, ci_str)
  } else {
    warning(paste("Expected BEDASSLE output file not found:", robj_file))
  }
}

res_df <- data.frame(Variable=vars, aE_aD_Mean=var_means, CI=var_CI)
out_csv <- paste0(args$out_prefix, "_BEDASSLE_RES_CI.csv")
write.csv(res_df, out_csv, row.names=FALSE)
cat(sprintf("Final Results CI table written to %s.\n", out_csv))
