#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
  library(dartR.popgen)
  library(parallel)
})

# ==========================================
# 1. Configuration and Parsing Module
# ==========================================
parse_cli_args <- function() {
  parser <- ArgumentParser(description="Create stairway plot blueprint files from a dartR VCF")
  
  # Required arguments
  parser$add_argument("-v", "--vcf", type="character", required=TRUE,
                      help="Path to input VCF file (can be gzipped)")
  parser$add_argument("-O", "--prefix", type="character", required=TRUE,
                      help="Prefix for output blueprint files")
  
  # Optional arguments - general
  parser$add_argument("-p", "--popmap", type="character", default=NULL,
                      help="Path to population file (CSV format). If skipped, simply makes the blueprint for the whole VCF instead")
  parser$add_argument("-o", "--outdir", type="character", default=".",
                      help="Output directory for blueprint files (default: .)")
  parser$add_argument("-s", "--stairway_plot_dir", type="character", default="stairway_plot_es",
                      help="Path to stairway_plot_es folder")
  
  # New Projection logic
  parser$add_argument("--projection", type="character", default="none",
                      help="Projection sampling strategy: 'none', 'conservative', 'max_sites', 'grid', or specific integer size. (default: 'none')")

  # Stairway plot customization
  parser$add_argument("--unfolded", action="store_false", dest="whether_folded",
                      help="SFS is unfolded (default is to use folded SFS)")
  parser$add_argument("--mu", type="double", default=1.2e-8,
                      help="Mutation rate per site per generation (default: 1.2e-8)")
  parser$add_argument("--read_length", type="integer", default=100,
                      help="Length of sequence reads (default: 100)")
  parser$add_argument("--year_per_generation", type="double", default=1.0,
                      help="Assumed generation time in years (default: 1)")
  parser$add_argument("--pct_training", type="double", default=0.67,
                      help="Percentage of sites for training (default: 0.67)")
  parser$add_argument("--nrand", type="character", default=NULL,
                      help="Comma-separated breakpoints for each try (e.g. '10,20,30'). Default is dynamically calculated based on sequence count.")
  parser$add_argument("--ninput", type="integer", default=200,
                      help="Number of input files to be created for each estimation (default: 200)")
  parser$add_argument("--xrange", type="character", default="0,0",
                      help="Time range formatted as 'xmin,xmax' (default: 0,0)")
  parser$add_argument("--yrange", type="character", default="0,0",
                      help="Ne range formatted as 'ymin,ymax' (default: 0,0)")
  parser$add_argument("--xspacing", type="double", default=2,
                      help="X axis spacing (default: 2)")
  parser$add_argument("--yspacing", type="double", default=2,
                      help="Y axis spacing (default: 2)")
  parser$add_argument("--fontsize", type="integer", default=12,
                      help="Font size (default: 12)")
  parser$add_argument("--smallest_size_of_sfs_bin", type="integer", default=1,
                      help="Smallest size of SFS bin (default: 1)")
  parser$add_argument("--largest_size_of_sfs_bin", type="integer", default=NULL,
                      help="Largest size of SFS bin")
  parser$add_argument("--random_seed", type="integer", default=NULL,
                      help="Random seed")
  
  args <- parser$parse_args()
  
  # Parse standard numeric vectors
  parse_numeric_vector <- function(x) {
    if (is.null(x)) return(NULL)
    as.numeric(unlist(strsplit(x, ",")))
  }
  
  args$nrand <- parse_numeric_vector(args$nrand)
  args$xrange <- parse_numeric_vector(args$xrange)
  args$yrange <- parse_numeric_vector(args$yrange)
  
  return(args)
}

# ==========================================
# 2. Mathematical & SFS Logic Module
# ==========================================

get_allele_counts <- function(gl) {
  gl_mat <- as.matrix(gl)
  k <- colSums(gl_mat, na.rm = TRUE)
  n_chr <- 2 * colSums(!is.na(gl_mat))
  list(k = k, n_chr = n_chr)
}

choose_projection <- function(n_chr, strategy = "none", original_nseq) {
  if (strategy == "none" || is.null(strategy)) {
    return("none")
  }

  if (!is.na(suppressWarnings(as.numeric(strategy)))) {
    return(as.numeric(strategy))
  }

  if (strategy == "conservative") {
    total_sites <- length(n_chr)
    candidates <- seq(2, original_nseq, by = 2)
    usable <- sapply(candidates, function(p) sum(n_chr >= p))
    # Strict 80% coverage
    valid_candidates <- candidates[usable >= 0.80 * total_sites]
    if (length(valid_candidates) > 0) {
      return(max(valid_candidates))
    } else {
      # Fallback if no projection retains 80%
      return(2)
    }
  }

  if (strategy == "max_sites") {
    candidates <- seq(2, original_nseq, by = 2)
    usable <- sapply(candidates, function(p) sum(n_chr >= p))
    score <- candidates * usable
    return(candidates[which.max(score)])
  }

  if (strategy == "grid") {
    return(seq(2, original_nseq, by = 2))
  }

  stop(sprintf("Unknown projection strategy: %s", strategy))
}

project_site <- function(k, n, m) {
  # using dhyper for robust fast combinations
  return(dhyper(0:m, k, n - k, m))
}

build_sfs <- function(k, n_chr, n_proj) {
  sfs <- rep(0, n_proj + 1)
  for (i in seq_along(k)) {
    if (!is.na(k[i]) && n_chr[i] >= n_proj) {
      sfs <- sfs + project_site(k[i], n_chr[i], n_proj)
    }
  }
  return(sfs)
}

fold_sfs <- function(sfs, n_proj) {
  # Expected sfs length is n_proj + 1. Index 1 = 0 mutant allele.
  folded <- sfs[1:(floor(n_proj / 2) + 1)]
  for (i in 2:length(folded)) {
    j <- n_proj - i + 2
    if (i != j) {
      folded[i] <- sfs[i] + sfs[j]
    }
  }
  return(folded)
}

compute_sfs <- function(dartR_vcf, args, original_nseq, counts, n_proj_val) {
  
  if (n_proj_val == "none") {
    # Wrapped standard dartR SFS calculation
    sfs_df <- gl.sfs(dartR_vcf, folded = args$whether_folded, singlepop = TRUE, minbinsize = args$smallest_size_of_sfs_bin, verbose = 0)
    sfs_vec <- c(data.frame(sfs_df)[, 1])
    
    usable_sites <- dartR_vcf@n.loc
    nseq_eff <- original_nseq
  } else {
    # Hypergeometric Downprojection Workflow
    unfolded_sfs <- build_sfs(counts$k, counts$n_chr, n_proj_val)
    if (args$whether_folded) {
      sfs_raw <- fold_sfs(unfolded_sfs, n_proj_val)
      sfs_vec <- round(sfs_raw[-1]) # Drop the k=0 bin (Stairway format starts with singletons)
    } else {
      sfs_vec <- round(unfolded_sfs[-c(1, length(unfolded_sfs))]) # Unfolded drops k=0 and k=n_proj
    }
    
    usable_sites <- sum(counts$n_chr >= n_proj_val)
    nseq_eff <- n_proj_val
  }
  
  # Max length safety warning
  max_sfs_length <- if (args$whether_folded) floor(nseq_eff / 2) else (nseq_eff - 1)
  if (length(sfs_vec) > max_sfs_length) {
    warning(sprintf("SFS length (%d) exceeds maximum expected (%d).", 
                    length(sfs_vec), max_sfs_length))
  }
  
  L_eff <- usable_sites * args$read_length
  
  return(list(sfs_vec = sfs_vec, nseq = nseq_eff, L = L_eff, usable_sites = usable_sites))
}

# ==========================================
# 3. File & Output Generation Module
# ==========================================
build_blueprint_content <- function(popid, nseq, L, sfs, nrand, args) {
  lines <- character(0)
  
  lines <- c(lines, "#Stairway Plot blueprint file")
  lines <- c(lines, "#input setting")
  lines <- c(lines, paste("popid:", popid, "# id of the population (no white space)"))
  lines <- c(lines, paste("nseq:", nseq, "# number of sequences"))
  lines <- c(lines, paste("L:", L, "# total number of observed nucleic sites, including polymorphic and monomorphic"))
  lines <- c(lines, paste("whether_folded:", tolower(args$whether_folded), "# whether the SFS is folded (true or false)"))
  
  sfs_str <- paste(sfs, collapse = "\t")
  lines <- c(lines, paste("SFS:", sfs_str, "# snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)"))
  
  if (!is.null(args$smallest_size_of_sfs_bin) && args$smallest_size_of_sfs_bin >= 1) {
    lines <- c(lines, paste("#smallest_size_of_SFS_bin_used_for_estimation:", args$smallest_size_of_sfs_bin, "# default is 1; to ignore singletons, uncomment this line and change this number to 2"))
  }
  
  if (!is.null(args$largest_size_of_sfs_bin) && args$largest_size_of_sfs_bin >= 1) {
    div_val <- ifelse(args$whether_folded, 2, 1)
    fold_str <- ifelse(args$whether_folded, "folded", "unfolded")
    lines <- c(lines, paste("#largest_size_of_SFS_bin_used_for_estimation:", args$largest_size_of_sfs_bin, sprintf("# default is nseq/%d for %s SFS", div_val, fold_str)))
  }
  
  lines <- c(lines, paste("pct_training:", args$pct_training, "# percentage of sites for training"))
  
  if (is.null(nrand)) {
    nrand <- c(round((nseq-2)/4, 0), round((nseq-2)/2, 0), round((nseq-2)*3/4, 0), round(nseq-2, 0))
  }
  nrand_str <- paste(nrand, collapse = "\t")
  lines <- c(lines, paste("nrand:", nrand_str, "# number of random break points for each try (separated by white space)"))
  
  lines <- c(lines, paste("project_dir:", file.path(args$outdir, popid), "# project directory"))
  lines <- c(lines, paste("stairway_plot_dir:", args$stairway_plot_dir, "# directory to the stairway plot files"))
  lines <- c(lines, paste("ninput:", args$ninput, "# number of input files to be created for each estimation"))
  
  if (!is.null(args$random_seed)) {
    lines <- c(lines, paste("#random_seed:", args$random_seed))
  }
  
  lines <- c(lines, "#output setting")
  lines <- c(lines, paste("mu:", args$mu, "# assumed mutation rate per site per generation"))
  lines <- c(lines, paste("year_per_generation:", args$year_per_generation, "# assumed generation time (in years)"))
  
  lines <- c(lines, "#plot setting")
  lines <- c(lines, paste("plot_title:", popid, "# title of the plot"))
  lines <- c(lines, paste("xrange:", paste(args$xrange, collapse = ","), "# Time (1k year) range; format: xmin,xmax; \"0,0\" for default"))
  lines <- c(lines, paste("yrange:", paste(args$yrange, collapse = ","), "# Ne (1k individual) range; format: ymin,ymax; \"0,0\" for default"))
  lines <- c(lines, paste("xspacing:", args$xspacing, "# X axis spacing"))
  lines <- c(lines, paste("yspacing:", args$yspacing, "# Y axis spacing"))
  lines <- c(lines, paste("fontsize:", args$fontsize, "# Font size"))
  
  return(lines)
}

write_blueprint <- function(output_file, content_lines) {
  writeLines(content_lines, output_file)
}

# ==========================================
# 4. Control Flow and Workflow Module
# ==========================================
generate_population_blueprint <- function(dartR_vcf, base_popid, args) {
  
  original_nseq <- length(dartR_vcf@ind.names) * 2
  counts <- get_allele_counts(dartR_vcf)
  
  n_proj_list <- choose_projection(counts$n_chr, strategy = args$projection, original_nseq = original_nseq)
  
  if (args$projection != "none") {
    plot_file <- file.path(args$outdir, paste0(base_popid, "_projection_plot.svg"))
    candidates <- seq(2, original_nseq, by = 2)
    usable <- sapply(candidates, function(p) sum(counts$n_chr >= p))
    
    df <- data.frame(Projected_Size = candidates, Sites_Retained = usable)
    
    p <- ggplot(df, aes(x = Projected_Size, y = Sites_Retained)) +
      geom_vline(xintercept = n_proj_list, color = "firebrick", linetype = "dashed", alpha = 0.8, linewidth = 0.8) +
      geom_line(color = "steelblue", linewidth = 1) +
      geom_point(color = "darkblue", size = 2) +
      theme_minimal() +
      theme(
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste("SFS Projection Trade-off:", base_popid),
        x = "Projected Sample Size (Chromosomes)",
        y = "Segregating Sites Retained"
      )
      
    suppressMessages(ggsave(plot_file, plot = p, width = 8, height = 6, device = "svg"))
    cat(sprintf("  --> Wrote projection plot: %s\n", plot_file))
  }
  
  for (n_proj_val in n_proj_list) {
    if (n_proj_val != "none" && n_proj_val < 2) {
        warning(sprintf("Calculated projection %s is too low, skipping...", n_proj_val))
        next
    }

    sfs_data <- compute_sfs(dartR_vcf, args, original_nseq, counts, n_proj_val)
    
    if (n_proj_val == "none") {
      proj_suffix <- "projNONE"
    } else {
      proj_suffix <- paste0("proj", n_proj_val)
    }
    
    popid_name <- paste0(base_popid, "_", proj_suffix)
    output_file <- file.path(args$outdir, paste0(popid_name, ".blueprint"))
    
    content <- build_blueprint_content(
      popid = popid_name, 
      nseq = sfs_data$nseq, 
      L = sfs_data$L, 
      sfs = sfs_data$sfs_vec, 
      nrand = args$nrand, 
      args = args
    )
    
    write_blueprint(output_file, content)
    cat(sprintf("  --> Wrote configuration: %s (Sites Used: %d / %d)\n", 
                output_file, sfs_data$usable_sites, length(counts$n_chr)))
  }
}

process_all_genotypes <- function(args) {
  
  cat("Reading VCF file:", args$vcf, "\n")
  dartR_vcf <- gl.read.vcf(args$vcf, mode='genotype')
  
  if (!is.null(args$popmap)) {
    
    cat("Reading population file:", args$popmap, "\n")
    pop_info <- read.csv(args$popmap)
    
    if (!"Indv" %in% colnames(pop_info)) {
      stop("Population file must contain an 'Indv' column with sample IDs matching VCF")
    }
    if (length(pop_info$Indv) != length(dartR_vcf@ind.names) || any(sort(pop_info$Indv) != sort(dartR_vcf@ind.names))) {
      stop("Population file sample names do not match VCF sample names")
    }
    
    hier_cols <- setdiff(colnames(pop_info), "Indv")
    if (length(hier_cols) == 0) {
      stop("Population file must contain at least one population column in addition to 'Indv'")
    }
    cat("Found population levels:", paste(hier_cols, collapse = ", "), "\n")
    
    row.names(pop_info) <- pop_info$Indv
    pop_info <- pop_info[dartR_vcf@ind.names, ]
    
    for (level in hier_cols) {
      cat("\nProcessing level:", level, "\n")
      pop(dartR_vcf) <- pop_info[, level]
      pops <- unique(pop_info[, level])
      
      for (pop_name in pops) {
        cat("Creating blueprint for population:", pop_name, "\n")
        
        # Subsetting to a single population by dropping all others
        singlepop_gl <- gl.drop.pop(dartR_vcf, pop.list = pops[pops != pop_name], mono.rm = TRUE, verbose = 0)
        nloc_orig <- singlepop_gl$n.loc
        singlepop_gl <- gl.filter.allna(singlepop_gl, verbose = 0)
        
        cat("Filtered", nloc_orig - singlepop_gl$n.loc, "loci with 100% missing data\n")
        
        # Configuration Details
        base_popid <- paste0(args$prefix, "_", level, "_", pop_name)
        generate_population_blueprint(singlepop_gl, base_popid, args)
      }
    }
  }
  
  # "All" Blueprint Generation
  cat("\nCreating blueprint for all samples combined\n")
  base_popid_all <- paste0(args$prefix, "_All")
  
  generate_population_blueprint(dartR_vcf, base_popid_all, args)
  
  cat("\nAll blueprint files created successfully!\n")
}

# ==========================================
# 5. Main Execution Block
# ==========================================
main <- function() {
  args <- parse_cli_args()
  
  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive = TRUE)
  }
  
  tryCatch({
    process_all_genotypes(args)
  }, error = function(e) {
    cat(sprintf("\nError occurred: %s\n", e$message))
    quit(status = 1)
  })
}

# Run the script when executed from terminal
if (!interactive()) {
  main()
}
