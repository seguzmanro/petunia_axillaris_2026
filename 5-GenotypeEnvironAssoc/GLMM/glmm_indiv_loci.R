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
})

pop_name_combns <- function(pop_names_vector){
  pop_name_combn_df <- data.frame(t(combn(pop_names_vector, 2)) ,stringsAsFactors = FALSE)
  names(pop_name_combn_df) <- c('Pop1', 'Pop2')
  return(pop_name_combn_df)
}

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

run_model <- function(dat_dataframe, list_of_vars, n_ger=2e6, burn=5e5, thin=750){
  formula <- as.formula(paste('gen', "~", paste(list_of_vars, collapse = " + ")))
  prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1)))
  model <- MCMCglmm(fixed = formula, random = ~ idv(mult.memb(~ Pop1 + Pop2)), data = dat_dataframe, prior = prior, pr=TRUE, nitt = n_ger, burnin = burn, thin = thin, verbose = FALSE)
  model
}

parser <- ArgumentParser(description='Fit individual loci MCMCglmm models')
parser$add_argument('--vcf', type='character', required=TRUE, help='Path to VCF file')
parser$add_argument('--popmap', type='character', required=TRUE, help='Population map CSV')
parser$add_argument('--env_dist', type='character', required=TRUE, help='CSV with pairwise geographic and environmental distances')
parser$add_argument('--out_prefix', type='character', required=TRUE, help='Prefix for output files')
parser$add_argument('--n_ger', type='integer', default=2000000, help='Number of MCMC generations')
parser$add_argument('--burnin', type='integer', default=500000, help='Number of burnin generations')
parser$add_argument('--thin', type='integer', default=750, help='Thinning interval')
parser$add_argument('--threads', type='integer', default=1, help='Number of threads')

args <- parser$parse_args()

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

# Extract CHROM_POS locus names from VCF for identifiable output
vcf_fix <- loaded_vcf@fix
locus_names <- paste(vcf_fix[, 'CHROM'], vcf_fix[, 'POS'], sep = '_')

outfreq_tables <- allele_frequency_tables(loaded_vcf, samples_info[,1], samples_info[,2])
outfreq_SO <- outfreq_tables[[1]]
outfreq_ST <- outfreq_tables[[2]]
names(outfreq_SO) <- locus_names
names(outfreq_ST) <- locus_names

model_names_simple <- c("1", env_var_names)
model_vars_simple <- as.list(c(1, env_var_names))
names(model_vars_simple) <- model_names_simple

cat(sprintf("Running individual loci models for %d loci with %d threads...\n", length(locus_names), args$threads))
cl <- makeCluster(args$threads, type = 'FORK')
registerDoParallel(cl)

DICs_all_list <- foreach(loc_name = locus_names) %dopar% {
  afd_dist <- AFD(outfreq_ST[[loc_name]])
  afd_dist <- afd_dist %>% rowwise() %>% mutate(
    p1 = min(as.character(Pop1), as.character(Pop2)),
    p2 = max(as.character(Pop1), as.character(Pop2))
  ) %>% ungroup()
  
  dat.loc <- merge(afd_dist, env_dist[, c("p1", "p2", env_var_names)], by=c("p1", "p2"))
  dat.loc$Pop1 <- factor(dat.loc$Pop1)
  dat.loc$Pop2 <- factor(dat.loc$Pop2)
  
  dat.loc$gen <- dat.loc$afd
  
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
