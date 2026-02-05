#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# Bayesian Outlier Detection Plotting Script
#
# A generalized, reusable script for analyzing and plotting results from both
# Bayescan and BayescEnv.
#
# Features:
# - MCMC convergence diagnostics
# - Outlier detection with FDR control
# - Visualization of results
# - Extraction of outlier loci from VCF files
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Required Libraries
#------------------------------------------------------------------------------#
required_packages <- c("coda", "readr", "vcfR", "ggmcmc")

# Install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}

# Load libraries
library(coda)
library(readr)
library(vcfR)
library(ggmcmc)

#------------------------------------------------------------------------------#
# Configuration
#------------------------------------------------------------------------------#
config <- list(
  # Analysis type: "bayescan" or "bayescenv"
  analysis_type = "bayescan",
  
  # Input paths
  working_dir = NULL,          # Auto-detect if NULL
  vcf_file = NULL,             # Path to VCF file
  results_dir = ".",           # Directory with result files
  sel_pattern = "\\.sel$",     # Pattern for .sel files
  fst_pattern = "_fst\\.txt$", # Pattern for _fst.txt files
  
  # Output settings
  output_prefix = NULL,        # Auto-detect from first result file
  fdr_threshold = 0.05,        # FDR threshold for outlier detection
  plot_height = 20,            # SVG plot height in inches
  plot_width = 20,             # SVG plot width in inches
  
  # MCMC settings
  thin_interval = 10,          # Thinning interval for MCMC chains
  
  # Plot settings
  point_size = 1,              # Size of points in scatter plot
  text_offset = 0.35,          # Offset for text labels
  highlight_loci = NULL,       # Loci to highlight (row names)
  name_highlighted = FALSE,    # Show names of highlighted loci
  add_text = TRUE              # Add text labels to outliers
)

#------------------------------------------------------------------------------#
# Helper Functions
#------------------------------------------------------------------------------#

#' Load and process MCMC chains from .sel files
#'
#' @param sel_files Vector of paths to .sel files
#' @param thin Thinning interval
#' @return List of mcmc objects
load_mcmc_chains <- function(sel_files, thin = 10) {
  chains <- list()
  
  for (i in seq_along(sel_files)) {
    chain <- read.table(sel_files[i], header = TRUE)
    # Remove first two columns (step and other metadata)
    chain <- chain[, -c(1, 2)]
    chains[[i]] <- mcmc(chain, thin = thin)
  }
  
  names(chains) <- sel_files
  return(chains)
}

#' Run MCMC convergence diagnostics
#'
#' @param chains List of mcmc objects
#' @param output_dir Directory to save diagnostics
run_convergence_diagnostics <- function(chains, output_dir = ".") {
  cat("Running MCMC convergence diagnostics...\n\n")
  
  # Process each chain separately
  for (i in seq_along(chains)) {
    chain_name <- names(chains)[i]
    base_name <- tools::file_path_sans_ext(basename(chain_name))
    
    cat(sprintf("Processing chain: %s\n", chain_name))
    
    # Plot traceplots
    plot(chains[[i]])
    
    # Save diagnostics using ggmcmc with chain-specific filename
    diag_file <- file.path(output_dir, paste0(base_name, "_diagnostics"))
    cat(sprintf("Saving detailed diagnostics to: %s.pdf\n", diag_file))
    ggmcmc(ggs(chains[[i]]), file = paste0(diag_file, ".pdf"), param_page = 5)
  }
  
  # Autocorrelation diagnosis
  cat("\nAutocorrelation diagnostics:\n")
  for (i in seq_along(chains)) {
    cat(sprintf("\nFile: %s\n", names(chains)[i]))
    print(autocorr.diag(chains[[i]]))
  }
  
  # Effective sample size
  cat("\nEffective sample sizes:\n")
  for (i in seq_along(chains)) {
    cat(sprintf("\nFile: %s\n", names(chains)[i]))
    print(effectiveSize(chains[[i]]))
  }
  
  # Geweke's convergence diagnostic
  cat("\nGeweke's convergence diagnostics:\n")
  for (i in seq_along(chains)) {
    cat(sprintf("\nFile: %s\n", names(chains)[i]))
    print(geweke.diag(chains[[i]], frac1 = 0.1, frac2 = 0.5))
  }
  
  # Heidelberg and Welch's convergence diagnostic
  cat("\nHeidelberg and Welch's convergence diagnostics:\n")
  for (i in seq_along(chains)) {
    cat(sprintf("\nFile: %s\n", names(chains)[i]))
    print(heidel.diag(chains[[i]], eps = 0.1, pvalue = 0.05))
  }
}

#' Plot results and detect outliers
#'
#' @param res Result table or file path
#' @param fdr FDR threshold
#' @param size Point size
#' @param pos Text offset
#' @param highlight Loci to highlight
#' @param name_highlighted Show names of highlighted loci
#' @param add_text Add text labels to outliers
#' @param analysis_type "bayescan" or "bayescenv"
#' @return List with outliers and count
plot_bayesian_results <- function(res, fdr = 0.05, size = 1, pos = 0.35, 
                                 highlight = NULL, name_highlighted = FALSE, 
                                 add_text = TRUE, analysis_type = "bayescan") {
  
  if (is.character(res)) {
    res <- read.table(res, header = TRUE)
  }
  
  # Determine column indices based on analysis type
  if (analysis_type == "bayescan") {
    col_fstat <- 5
  } else if (analysis_type == "bayescenv") {
    col_fstat <- 7
  } else {
    stop(sprintf("Unknown analysis type: %s. Use 'bayescan' or 'bayescenv'.", analysis_type))
  }
  
  col_q <- col_fstat - 2
  
  # Identify highlight and non-highlight rows
  highlight_rows <- if (!is.null(highlight)) {
    which(is.element(as.numeric(rownames(res)), highlight))
  } else {
    integer(0)
  }
  non_highlight_rows <- setdiff(1:nrow(res), highlight_rows)
  
  # Detect outliers
  outliers <- as.integer(rownames(res[res[, col_q] <= fdr, ]))
  ok_outliers <- length(outliers) > 0
  
  # Cap very small q-values for plotting
  res[res[, col_q] <= 0.0001, col_q] <- 0.0001
  
  # Create plot
  plot(log10(res[, col_q]), res[, col_fstat],
       xlim = rev(range(log10(res[, col_q]))),
       xlab = "log10(q value)",
       ylab = names(res[col_fstat]),
       type = "n")
  points(log10(res[non_highlight_rows, col_q]), res[non_highlight_rows, col_fstat],
         pch = 19, cex = size)
  
  if (name_highlighted) {
    if (length(highlight_rows) > 0) {
      text(log10(res[highlight_rows, col_q]), res[highlight_rows, col_fstat],
           rownames(res[highlight_rows, ]), col = "red", cex = size * 1.2, font = 2)
    }
  } else {
    if (length(highlight_rows) > 0) {
      points(log10(res[highlight_rows, col_q]), res[highlight_rows, col_fstat],
             col = "red", pch = 19, cex = size)
    }
    if (ok_outliers && add_text) {
      text(log10(res[res[, col_q] <= fdr, ][, col_q]) + 
             pos * (round(runif(nrow(res[res[, col_q] <= fdr, ]), 1, 2)) * 2 - 3),
           res[res[, col_q] <= fdr, ][, col_fstat],
           rownames(res[res[, col_q] <= fdr, ]), cex = size)
    }
  }
  lines(c(log10(fdr), log10(fdr)), c(-1, 1), lwd = 2)
  
  return(list("outliers" = outliers, "nb_outliers" = length(outliers)))
}

#' Extract outlier loci from VCF file
#'
#' @param vcf_path Path to VCF file
#' @param outliers List of outlier indices (row numbers)
#' @return Data frame with CHROM and POS columns for outliers
extract_outliers_from_vcf <- function(vcf_path, outliers) {
  if (is.null(vcf_path)) {
    warning("VCF file path not specified. Skipping outlier extraction.")
    return(NULL)
  }
  
  if (!file.exists(vcf_path)) {
    stop(sprintf("VCF file not found: %s", vcf_path))
  }
  
  vcf <- read.vcfR(vcf_path)
  return(vcf@fix[outliers, c(1, 2)])
}

#' Process all result files
#'
#' @param config Configuration list
process_results <- function(config) {
  # Set working directory
  if (!is.null(config$working_dir)) {
    setwd(config$working_dir)
  }
  
  cat(sprintf("Working directory: %s\n", getwd()))
  
  # Find result files
  sel_files <- list.files(path = config$results_dir, 
                          pattern = config$sel_pattern, 
                          full.names = TRUE)
  fst_files <- list.files(path = config$results_dir, 
                          pattern = config$fst_pattern, 
                          full.names = TRUE)
  
  cat(sprintf("Found %d .sel files and %d _fst.txt files\n", 
              length(sel_files), length(fst_files)))
  
  # Load and process MCMC chains
  if (length(sel_files) > 0) {
    chains <- load_mcmc_chains(sel_files, config$thin_interval)
    run_convergence_diagnostics(chains, output_dir = dirname(sel_files[1]))
  }
  
  # Process Fst result files
  if (length(fst_files) > 0) {
    # Generate output prefixes
    res_files <- lapply(fst_files, function(x) {
      sub(config$fst_pattern, "", x)
    })
    
    # Load VCF file if specified
    if (!is.null(config$vcf_file)) {
      if (!file.exists(config$vcf_file)) {
        warning(sprintf("VCF file not found: %s. Will skip outlier extraction.", config$vcf_file))
        vcf <- NULL
      } else {
        cat(sprintf("Loading VCF file: %s\n", config$vcf_file))
        vcf <- read.vcfR(config$vcf_file)
      }
    } else {
      vcf <- NULL
    }
    
    # Process each result file
    for (i in seq_along(fst_files)) {
      cat(sprintf("\nProcessing file: %s\n", fst_files[i]))
      
      # Create output filename prefix
      output_prefix <- if (!is.null(config$output_prefix)) {
        config$output_prefix
      } else {
        basename(res_files[[i]])
      }
      
      # Generate plot filename
      plot_file <- sprintf("%s_Outliers_%s_%02d.svg", 
                          res_files[[i]],
                          toupper(substr(config$analysis_type, 1, 2)),
                          config$fdr_threshold * 100)
      
      # Generate output table filename
      table_file <- sprintf("%s_Outliers_%s_%02d.txt", 
                           res_files[[i]],
                           toupper(substr(config$analysis_type, 1, 2)),
                           config$fdr_threshold * 100)
      
      # Create plot
      cat("Creating plot...\n")
      svg(plot_file, height = config$plot_height, width = config$plot_width)
      outliers <- plot_bayesian_results(
        fst_files[i],
        fdr = config$fdr_threshold,
        size = config$point_size,
        pos = config$text_offset,
        highlight = config$highlight_loci,
        name_highlighted = config$name_highlighted,
        add_text = config$add_text,
        analysis_type = config$analysis_type
      )
      dev.off()
      
      cat(sprintf("Outliers detected: %d\n", outliers$nb_outliers))
      cat(sprintf("Plot saved to: %s\n", plot_file))
      
      # Save outliers to file
      if (!is.null(vcf) && outliers$nb_outliers > 0) {
        cat("Extracting outlier loci from VCF file...\n")
        outlier_loci <- vcf@fix[outliers$outliers, c(1, 2)]
        write.table(outlier_loci, file = table_file, 
                    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
        cat(sprintf("Outliers table saved to: %s\n", table_file))
      }
    }
  }
}

#' Print usage information
print_usage <- function() {
  cat(
    "Bayesian Outlier Detection Plotting Script\n\n",
    "Usage:\n",
    "1. Set working directory and VCF file path in the script\n",
    "2. Run the script with `Rscript plot_bayes.R`\n",
    "3. Results will be saved in the specified working directory\n\n",
    "For Bayescan:\n",
    "- .sel files should contain MCMC chain information\n",
    "- _fst.txt files should contain result data\n\n",
    "For BayescEnv:\n",
    "- Set analysis_type to 'bayescenv'\n",
    "- Adjust patterns if filenames have different formats\n\n"
  )
}

#------------------------------------------------------------------------------#
# Main Execution
#------------------------------------------------------------------------------#

main <- function() {
  cat("Bayesian Outlier Detection Plotting Script\n")
  cat("==========================================\n\n")
  
  # Example configurations
  
  # Uncomment this for Bayescan
  # config$analysis_type <- "bayescan"
  # config$working_dir <- "./"
  # config$vcf_file <- "../../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf"
  # config$results_dir <- "."
  # config$fst_pattern <- "_fst\\.txt$"
  # config$fdr_threshold <- 0.01
  
  # Uncomment this for BayescEnv (adjust paths as needed)
  # config$analysis_type <- "bayescenv"
  # config$working_dir <- "../../5-GenotypeEnvironAssoc/Bayescenv/"
  # config$vcf_file <- "../../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf.gz"
  # config$results_dir <- "../../5-GenotypeEnvironAssoc/Bayescenv/Results"
  # config$fst_pattern <- "_fst\\.txt$"  # Adjust if filenames are different
  # config$fdr_threshold <- 0.01
  
  # Run processing
  tryCatch({
    process_results(config)
  }, error = function(e) {
    cat(sprintf("\nError: %s\n", e$message))
    cat("\nStack trace:\n")
    print(e)
    cat(sprintf("\nCurrent working dir: %s\n", getwd()))
    cat(sprintf("\nCurrent config params: %s\n", config))
    cat("\n")
    print_usage()
  })
  
  cat("\nProcessing completed.\n")
}

# Run main if script is executed directly
if (sys.nframe() == 0) {
  main()
}

#------------------------------------------------------------------------------#
# End of Script
#------------------------------------------------------------------------------#
