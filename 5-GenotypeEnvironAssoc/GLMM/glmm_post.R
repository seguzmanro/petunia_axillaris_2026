#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MCMCglmm)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(argparse)
})

## Functions

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

post_extract_fixed_effects <- function(model){
  # model$Sol includes random effect samples when pr=TRUE; restrict to fixed effects only
  fixed_cols <- colnames(model$X)
  sol <- model$Sol[, fixed_cols, drop=FALSE]
  
  res <- data.frame(
    term = colnames(sol),
    post_mean = apply(sol, 2, mean),
    lower_95 = apply(sol, 2, quantile, 0.025),
    upper_95 = apply(sol, 2, quantile, 0.975),
    pMCMC = 2 * pmin(
      colMeans(sol > 0),
      colMeans(sol < 0)
    )
  )
  
  return(res)
}

post_compute_variance_partition <- function(model){
  X <- model$X
  # model$Sol includes random effect samples when pr=TRUE; restrict to fixed effects only
  beta <- model$Sol[, colnames(X), drop=FALSE]
  
  fitted_vals <- X %*% t(beta)
  varF <- mean(apply(fitted_vals, 2, var))
  
  vcv_means <- colMeans(model$VCV)
  
  var_random <- sum(vcv_means[grep("Pop", names(vcv_means))])
  var_residual <- sum(vcv_means[grep("units", names(vcv_means))])
  
  total_var <- varF + var_random + var_residual
  
  R2_marginal <- varF / total_var
  R2_conditional <- (varF + var_random) / total_var
  
  return(list(
    var_fixed = varF,
    var_random = var_random,
    var_residual = var_residual,
    R2_marginal = R2_marginal,
    R2_conditional = R2_conditional
  ))
}

parser <- ArgumentParser(description='Generate summary CSV from fitted MCMCglmm models')
parser$add_argument('--models_file', type='character', required=TRUE, help='Path to _fst_models.RData file')
parser$add_argument('--dics_file', type='character', required=FALSE, default=NULL, help='Path to _fst_DICs_res.csv (optional, DIC computed from models if not provided)')
parser$add_argument('--out_prefix', type='character', required=TRUE, help='Prefix for output CSV')

args <- parser$parse_args()

load(args$models_file)

DICs_table <- make_DICs_result_table(models_rad, names(models_rad))

all_fixed <- list()
all_var <- list()

for(model_name in names(models_rad)){
  model <- models_rad[[model_name]]
  
  fixed_df <- post_extract_fixed_effects(model)
  fixed_df$model <- model_name
  
  var_stats <- post_compute_variance_partition(model)
  
  var_df <- data.frame(
    model = model_name,
    var_fixed = var_stats$var_fixed,
    var_random = var_stats$var_random,
    var_residual = var_stats$var_residual,
    R2_marginal = var_stats$R2_marginal,
    R2_conditional = var_stats$R2_conditional
  )
  
  all_fixed[[model_name]] <- fixed_df
  all_var[[model_name]] <- var_df
}

fixed_results <- bind_rows(all_fixed)
var_results <- bind_rows(all_var)

fixed_wide <- fixed_results %>%
  pivot_wider(
    id_cols = model,
    names_from = term,
    values_from = c(post_mean, lower_95, upper_95, pMCMC),
    names_sep = "_"
  )

combined <- var_results %>%
  left_join(DICs_table %>% tibble::rownames_to_column("model"), by = "model") %>%
  left_join(fixed_wide, by = "model")

write.csv(combined, paste0(args$out_prefix, "_glmm_full_results.csv"), row.names=FALSE)

cat("Post-analysis summary saved.\n")
