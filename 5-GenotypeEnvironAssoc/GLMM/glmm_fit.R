#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(hierfstat)
  library(MCMCglmm)
  library(introgress)
  library(reshape2)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(argparse)
  library(coda)
  library(ggmcmc)
})

## Functions

pop_name_combns <- function(pop_names_vector){
  pop_name_combn_df <- data.frame(t(combn(pop_names_vector, 2)) ,stringsAsFactors = FALSE)
  names(pop_name_combn_df) <- c('Pop1', 'Pop2')
  return(pop_name_combn_df)
}

make_model_names <- function(model_var_names){
  vars_combn_2 <- combn(model_var_names, 2)
  model_vars <- list(1)

  for (i in (1:length(model_var_names))){
    model_vars[[length(model_vars)+1]] <- model_var_names[i]
  }

  for (i in (1:ncol(vars_combn_2))){
    if ('geog' %in% vars_combn_2[,i]){
      model_vars[[length(model_vars)+1]] <- vars_combn_2[,i]
    }
  }

  model_names <- list()
  for (i in (1:length(model_vars))){
    model_names[[i]] <- paste0(model_vars[[i]], collapse='_')
  }
  return(list(model_vars, model_names))
}

run_model <- function(dat_dataframe, list_of_vars, n_ger=2e6, burn=5e5, thin=750){
  formula <- as.formula(paste('gen', "~", paste(list_of_vars, collapse = " + ")))
  model <- MCMCglmm(fixed = formula, random = ~ idv(mult.memb(~ Pop1 + Pop2)), data = dat_dataframe, pr=TRUE, nitt = n_ger, burnin = burn, thin = thin, verbose = FALSE)
  model
}

run_diagnostics <- function(model_list, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ess_results <- list()
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    if (!inherits(model, "MCMCglmm")) {
      warning(paste("Object", model_name, "is not an MCMCglmm model. Skipping."))
      next
    }

    mcmc_chain <- as.mcmc(cbind(model$Sol, model$VCV))

    ess_values <- data.frame(effectiveSize(mcmc_chain))
    ess_file <- file.path(output_dir, paste0(model_name, "_ess.csv"))
    write.csv(ess_values, ess_file)
    cat("Saved ESS table for model:", model_name, "to:", ess_file, "\n")
    ess_results[[model_name]] <- ess_values

    pdf_file <- file.path(output_dir, paste0(model_name, "_diagnostics.pdf"))
    ggmcmc(ggs(mcmc_chain), file = pdf_file, param_page = 5)
    cat("Saved diagnostic plots for model:", model_name, "to:", pdf_file, "\n")
  }
  return(ess_results)
}

parser <- ArgumentParser(description='Fit MCMCglmm models for genotype-environment association')
parser$add_argument('--vcf', type='character', required=TRUE, help='Path to VCF file')
parser$add_argument('--popmap', type='character', required=TRUE, help='Population map CSV. First col=sample names, second col=population')
parser$add_argument('--env_dist', type='character', required=TRUE, help='CSV with pairwise geographic and environmental distances')
parser$add_argument('--out_prefix', type='character', required=TRUE, help='Prefix for output files')
parser$add_argument('--n_ger', type='integer', default=2000000, help='Number of MCMC generations')
parser$add_argument('--burnin', type='integer', default=500000, help='Number of burnin generations')
parser$add_argument('--thin', type='integer', default=750, help='Thinning interval')
parser$add_argument('--threads', type='integer', default=1, help='Number of threads')

args <- parser$parse_args()

### load data
loaded_vcf <- read.vcfR(args$vcf)
samples_info <- read.csv(args$popmap)
rownames(samples_info) <- samples_info[,1]

env_dist <- read.csv(args$env_dist)
env_dist <- env_dist %>% rowwise() %>% mutate(
  p1 = min(as.character(pop1), as.character(pop2)),
  p2 = max(as.character(pop1), as.character(pop2))
) %>% ungroup()
env_var_names <- setdiff(colnames(env_dist), c("pop1", "pop2", "p1", "p2"))

out_prefix <- args$out_prefix

cat("Computing Pairwise WC Fst...\n")
rad_dist <- pairwise.WCfst(vcfR2genind(loaded_vcf, pop=samples_info[,2]))
rad_dist <- rad_dist[unique(samples_info[,2]), unique(samples_info[,2])]

gen_dat <- pop_name_combns(rownames(rad_dist))
gen_dat$gen <- apply(gen_dat, 1, function(x){rad_dist[x[1], x[2]]})
gen_dat$gen[gen_dat$gen < 0] <- 0
gen_dat$gen <- scale(gen_dat$gen)
gen_dat <- gen_dat %>% rowwise() %>% mutate(
  p1 = min(as.character(Pop1), as.character(Pop2)),
  p2 = max(as.character(Pop1), as.character(Pop2))
) %>% ungroup()

dat.rad <- merge(gen_dat, env_dist[, c("p1", "p2", env_var_names)], by=c("p1", "p2"))
dat.rad$Pop1 <- factor(dat.rad$Pop1)
dat.rad$Pop2 <- factor(dat.rad$Pop2)

model_objs <- make_model_names(env_var_names)
model_vars <- model_objs[[1]]
model_names <- model_objs[[2]]

cat(sprintf("Running %d Fst models with %d threads...\n", length(model_vars), args$threads))
cl <- makeCluster(args$threads, type = 'FORK')
registerDoParallel(cl)

models_rad <- foreach(iteration=(1:length(model_vars))) %dopar% {
  run_model(dat.rad, model_vars[[iteration]], n_ger = args$n_ger, burn = args$burnin, thin = args$thin)
}
stopCluster(cl)

names(models_rad) <- model_names

save(models_rad, file=paste0(out_prefix, '_fst_models.RData'))
cat("Fst modeling complete.\n")

cat("Running MCMC diagnostics for FST models...\n")
run_diagnostics(models_rad, paste0(out_prefix, "_mcmc_plots"))
