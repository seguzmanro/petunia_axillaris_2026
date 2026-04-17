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

# allele frequency differential function
AFD <- function(allfreqmat){
  afdmat <- pop_name_combns(colnames(allfreqmat))
  afdmat$afd <- rep(NA, nrow(afdmat))
  for(i in 1 : nrow(afdmat)){
    afdmat[i,'afd'] <- delta(allfreqmat[,afdmat[i,1]],allfreqmat[,afdmat[i,2]])
  }
  return(afdmat)
}

allele_frequency_tables <- function(vcf, sample_order=NULL, population_vector){
  gen <- vcfR2genind(vcf)
  if (!is.null(sample_order)){
    gen@tab <- gen@tab[sample_order,]
  }
  genp <- genind2genpop(gen, population_vector)
  allele_count_pop <- genp@tab
  odd_alleles <- seq(1, ncol(allele_count_pop),2)
  even_alleles <- seq(2, ncol(allele_count_pop),2)
  outfreq_SO <- list()
  outfreq_ST <- list()
  for (i in 1:length(even_alleles)){
    outfreq_SO[[i]] <- t(allele_count_pop[,c(odd_alleles[i],even_alleles[i])])
    outfreq_ST[[i]] <- matrix(nrow = nrow(outfreq_SO[[i]]), ncol = ncol(outfreq_SO[[i]]))
    outfreq_ST[[i]][1,] <- outfreq_SO[[i]][1,]/(outfreq_SO[[i]][1,]+outfreq_SO[[i]][2,])
    outfreq_ST[[i]][2,] <- outfreq_SO[[i]][2,]/(outfreq_SO[[i]][1,]+outfreq_SO[[i]][2,])
    colnames(outfreq_ST[[i]]) <- colnames(outfreq_SO[[i]])
  }
  return(list(outfreq_SO, outfreq_ST))
}

make_model_names <- function(model_var_names){
  vars_combn_2 <- combn(model_var_names, 2) 
  model_vars <- list(1)
  
  # all the simple 1-variable models
  for (i in (1:length(model_var_names))){
    model_vars[[length(model_vars)+1]] <- model_var_names[i]
  }
  
  # 2-variable models including geog
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

make_DICs_result_table <- function(list_of_models, names_of_target_models){
  DICs_table <- data.frame(array(NA, c(length(names_of_target_models), 3)))
  row.names(DICs_table) <- names_of_target_models
  colnames(DICs_table) <- c('DIC', 'deltaDIC', 'DICweight')
  
  for (i in (1:length(names_of_target_models))){
    for (j in (1:length(list_of_models))){
      if (names(list_of_models)[j] == row.names(DICs_table)[[i]]){
        DICs_table[i,'DIC'] <- list_of_models[[j]]$DIC
      }
    }
  }
  
  DICs_table[,'deltaDIC'] <- with(data.frame(DICs_table), DIC-min(DIC))
  DICs_table[,'DICweight'] <- with(data.frame(DICs_table), exp(-deltaDIC / 2) / sum(exp(-deltaDIC / 2)))
  return(DICs_table)
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
    
    # Calculate ESS
    ess_values <- data.frame(effectiveSize(mcmc_chain))
    ess_file <- file.path(output_dir, paste0(model_name, "_ess.csv"))
    write.csv(ess_values, ess_file)
    cat("Saved ESS table for model:", model_name, "to:", ess_file, "\n")
    ess_results[[model_name]] <- ess_values
    
    # Plot diagnostics
    pdf_file <- file.path(output_dir, paste0(model_name, "_diagnostics.pdf"))
    ggmcmc(ggs(mcmc_chain), file = pdf_file, param_page = 5)
    cat("Saved diagnostic plots for model:", model_name, "to:", pdf_file, "\n")
  }
  return(ess_results)
}

parser <- ArgumentParser(description='Refactored GLMM R script for Pop Genomics')
parser$add_argument('--vcf', type='character', required=TRUE, help='Path to VCF file')
parser$add_argument('--popmap', type='character', required=TRUE, help='Population map CSV. First col=sample names, second col=population')
parser$add_argument('--env_dist', type='character', required=TRUE, help='CSV file with pairwise geographic and environmental distances (e.g. pop1, pop2, geog, env_var1...)')
parser$add_argument('--indiv_loci', action='store_true', help='If selected, will run Fst model AND additionally run AFD model for individual loci.')
parser$add_argument('--out_prefix', type='character', required=TRUE, help='Prefix for output files (can include dir path)')
parser$add_argument('--n_ger', type='integer', default=2000000, help='Number of MCMC generations, default=2e6')
parser$add_argument('--burnin', type='integer', default=500000, help='Number of burnin generations, default=5e5')
parser$add_argument('--thin', type='integer', default=750, help='Number of generations thinned, default=750')
parser$add_argument('--threads', type='integer', default=1, help='Number of threads')

args <- parser$parse_args()

### load data
loaded_vcf <- read.vcfR(args$vcf)
samples_info <- read.csv(args$popmap)
rownames(samples_info) <- samples_info[,1]

# load environment/geographic pairwise distance
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

# Get models
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
DICs.rad_full <- make_DICs_result_table(models_rad, names(models_rad))

save(models_rad, file=paste0(out_prefix, '_fst_models.RData'))
write.csv(DICs.rad_full, paste0(out_prefix, '_fst_DICs_res.csv'), row.names=TRUE)
cat("Fst modeling complete.\n")

cat("Running MCMC diagnostics for FST models...\n")
run_diagnostics(models_rad, paste0(out_prefix, "_mcmc_plots"))

if (args$indiv_loci) {
  cat("Running individual loci models...\n")
  outfreq_tables <- allele_frequency_tables(loaded_vcf, samples_info[,1], samples_info[,2])
  outfreq_SO <- outfreq_tables[[1]]
  outfreq_ST <- outfreq_tables[[2]]
  names(outfreq_SO) <- c(paste('loc', 1:length(outfreq_SO) ,sep = '_'))
  names(outfreq_ST) <- c(paste('loc', 1:length(outfreq_ST) ,sep = '_'))
  
  outfrequniqs <- lapply(outfreq_SO, function(x) {
    tmp <- x
    tmp[tmp[,1] / tmp[,2] == 0, 2] <- 1
    tmp[tmp[,2] / tmp[,1] == 0, 1] <- 1
    tmp
  })
  Ndivpops <- lapply(outfrequniqs, function(x) { nrow(unique(x)) })
  cutoff <- 1
  loci <- names(Ndivpops)[Ndivpops >= cutoff]
  
  # For individual loci we only test the null model intercept + 1-var models
  model_names_simple <- c("1", env_var_names)
  model_vars_simple <- as.list(c(1, env_var_names))
  names(model_vars_simple) <- model_names_simple
  
  cl <- makeCluster(args$threads, type = 'FORK')
  registerDoParallel(cl)

  DICs_all_list <- foreach(iteration=(1:length(loci))) %dopar% {
    loc_name <- loci[iteration]
    afd_dist <- AFD(outfreq_ST[[loc_name]])
    afd_dist <- afd_dist %>% rowwise() %>% mutate(
      p1 = min(as.character(Pop1), as.character(Pop2)),
      p2 = max(as.character(Pop1), as.character(Pop2))
    ) %>% ungroup()
    
    dat.loc <- merge(afd_dist, env_dist[, c("p1", "p2", env_var_names)], by=c("p1", "p2"))
    dat.loc$Pop1 <- factor(dat.loc$Pop1)
    dat.loc$Pop2 <- factor(dat.loc$Pop2)
    
    # Needs a gen column
    dat.loc$gen <- dat.loc$afd
    # afd wasn't originally scaled in the script
    
    loc_DICs <- data.frame(model=model_names_simple, DIC=NA, deltaDIC=NA, DICweight=NA)
    
    for (v_idx in seq_along(model_names_simple)){
      mod <- run_model(dat.loc, list_of_vars=model_vars_simple[[v_idx]], n_ger = args$n_ger, burn = args$burnin, thin=args$thin)
      loc_DICs$DIC[v_idx] <- mod$DIC
    }
    
    loc_DICs$deltaDIC <- loc_DICs$DIC - min(loc_DICs$DIC)
    loc_DICs$DICweight <- exp(-loc_DICs$deltaDIC / 2) / sum(exp(-loc_DICs$deltaDIC / 2))
    
    list(loc=loc_name, results=loc_DICs)
  }
  stopCluster(cl)
  
  save(DICs_all_list, file=paste0(out_prefix, '_indiv_loci_DICs.RData'))
  cat("Individual loci modeling complete.\n")
}
