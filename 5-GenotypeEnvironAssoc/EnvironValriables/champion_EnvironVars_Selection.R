library(argparse)
library(dplyr)
library(tidyr)
library(purrr)
library(terra)

# Argument parsing
parser <- ArgumentParser(description='Champion strategy for Environmental variable selection')
parser$add_argument('--popmap', type="character", required=TRUE, help='Path to the population map CSV')
parser$add_argument('--env_table', type="character", required=TRUE, help='Path to the environmental table CSV')
parser$add_argument('--output', type="character", required=TRUE, help='Output file name for the selected variables distances')
parser$add_argument('--out_dir_txt', type="character", default=NULL, help='Optional directory to write Bayescenv txt files')
parser$add_argument('--dendrogram_plot', type="character", required=TRUE, help='Path to output dendrogram plot')
parser$add_argument('--rankings_out', type="character", required=TRUE, help='Path to output variable PCA rankings')
parser$add_argument('--num_clusters_soil', type="integer", default=2, help='Number of soil clusters')
parser$add_argument('--num_clusters_atmos', type="integer", default=3, help='Number of atmos clusters')
parser$add_argument('--num_champions_per_cluster', type="integer", default=1, help='Number of champions to select per cluster automatically')
parser$add_argument('--champions', type="character", nargs='+', default="auto", help='List of champion variables or "auto"')

args <- parser$parse_args()

# Helper function for pairwise distances
pairwise_df <- function(mat) {
  as.data.frame(as.table(mat)) %>%
    setNames(c("pop1", "pop2", "distance")) %>%
    mutate(
      pop1 = as.character(pop1),
      pop2 = as.character(pop2),
      pop_lo = pmin(pop1, pop2),
      pop_hi = pmax(pop1, pop2)
    ) %>%
    filter(pop_lo != pop_hi) %>%
    distinct(pop_lo, pop_hi, .keep_all = TRUE) %>%
    transmute(
      pop1 = pop_lo,
      pop2 = pop_hi,
      distance
    )
}

# 1. Read data
eco_var_table <- read.csv(args$env_table, stringsAsFactors = TRUE)
paxil_popmap <- read.csv(args$popmap)
rownames(eco_var_table) <- eco_var_table$Pop

# Keep only populations in popmap
eco_var_table <- eco_var_table[paxil_popmap[,'Pop'][!duplicated(paxil_popmap[,'Pop'])],]

# Geographical distances
geog <- vect(
  eco_var_table[, c("lon", "lat")],
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

# Extract environ only
environ_only <- eco_var_table[, grepl("mean|CM|wc|Alt", names(eco_var_table))]

# Separate soil and atmos variables
soil_vars <- names(environ_only)[grepl("mean", names(environ_only))]
atmos_vars <- names(environ_only)[grepl("CM|wc|Alt", names(environ_only))]

process_cluster_group <- function(vars, num_clusters, group_name) {
  if (length(vars) == 0) return(list(champions = c(), rankings = list(), dendro = NULL, num_clusters = num_clusters))
  
  group_data <- environ_only[, vars, drop=FALSE]
  cor_matrix <- abs(cor(group_data, use = "pairwise.complete.obs", method = "pearson"))
  dist_matrix <- as.dist(1 - cor_matrix)
  env_cluster <- hclust(dist_matrix, method = "ward.D2")
  
  cluster_assignments <- cutree(env_cluster, k = num_clusters)
  cluster_list <- split(names(cluster_assignments), cluster_assignments)
  
  local_champion_vars <- c()
  rankings_list <- list()
  is_auto <- (length(args$champions) == 1 && args$champions[1] == "auto")
  
  for (i in seq_along(cluster_list)) {
    cluster_vars <- cluster_list[[i]]
    cluster_data <- environ_only[, cluster_vars, drop=FALSE]
    cluster_data <- na.omit(cluster_data)
    
    if (length(cluster_vars) > 1) {
      var_check <- apply(cluster_data, 2, var)
      valid_vars <- names(var_check[var_check > 0])
      
      if (length(valid_vars) > 1) {
        pca <- prcomp(cluster_data[, valid_vars], scale. = TRUE, center = TRUE)
        loadings_pc1 <- abs(pca$rotation[, 1])
        ordered_vars <- sort(loadings_pc1, decreasing = TRUE)
        
        constant_vars <- setdiff(cluster_vars, valid_vars)
        if (length(constant_vars) > 0) {
          zero_loadings <- rep(0, length(constant_vars))
          names(zero_loadings) <- constant_vars
          ordered_vars <- c(ordered_vars, zero_loadings)
        }
      } else {
        ordered_vars <- rep(0, length(cluster_vars))
        names(ordered_vars) <- cluster_vars
        if (length(valid_vars) == 1) ordered_vars[valid_vars[1]] <- 1.0
      }
      
      if (is_auto) {
        num_to_pick <- min(args$num_champions_per_cluster, length(ordered_vars))
        champ <- names(ordered_vars)[1:num_to_pick]
        local_champion_vars <- c(local_champion_vars, champ)
      }
      
      rankings_list[[i]] <- data.frame(
        Group = group_name,
        Cluster = i,
        Variable = names(ordered_vars),
        PC1_Absolute_Loading = ordered_vars,
        row.names = NULL
      )
      
    } else {
      if (is_auto) {
        champ <- cluster_vars[1]
        local_champion_vars <- c(local_champion_vars, champ)
      }
      rankings_list[[i]] <- data.frame(
        Group = group_name,
        Cluster = i,
        Variable = cluster_vars[1],
        PC1_Absolute_Loading = 1,
        row.names = NULL
      )
    }
  }
  
  return(list(
    champions = local_champion_vars,
    rankings_list = rankings_list,
    dendro = env_cluster,
    num_clusters = num_clusters
  ))
}

is_auto <- (length(args$champions) == 1 && args$champions[1] == "auto")
if (is_auto) cat("Running automated PCA-based champion variable selection...\n") else cat("Using manually provided champion variables...\n")

cat("Processing Soil variables...\n")
soil_res <- process_cluster_group(soil_vars, args$num_clusters_soil, "Soil")

cat("Processing Atmospheric/Altitudinal variables...\n")
atmos_res <- process_cluster_group(atmos_vars, args$num_clusters_atmos, "Atmos_Alt")

# Merge results
pdf(args$dendrogram_plot, width = 10, height = 7)
if (!is.null(soil_res$dendro)) {
  plot(soil_res$dendro, main = "Hierarchical Clustering of Soil Variables", xlab = "Variables", ylab = "Distance (1 - |r|)", sub = "", cex = 0.8)
  rect.hclust(soil_res$dendro, k = soil_res$num_clusters, border = "red")
}
if (!is.null(atmos_res$dendro)) {
  plot(atmos_res$dendro, main = "Hierarchical Clustering of Atmos/Alt Variables", xlab = "Variables", ylab = "Distance (1 - |r|)", sub = "", cex = 0.8)
  rect.hclust(atmos_res$dendro, k = atmos_res$num_clusters, border = "red")
}
invisible(dev.off())

if (is_auto) {
  champion_vars <- c(soil_res$champions, atmos_res$champions)
} else {
  champion_vars <- args$champions
}

# Attach Is_Champion to rankings
all_rankings_list <- c(soil_res$rankings_list, atmos_res$rankings_list)
all_rankings <- bind_rows(all_rankings_list) %>%
  mutate(Is_Champion = (Variable %in% champion_vars))

write.csv(all_rankings, args$rankings_out, row.names = FALSE)
cat("Wrote cluster rankings to", args$rankings_out, "\n")
cat("Selected champions:\n", paste(champion_vars, collapse = ", "), "\n")

# Validate user champions exist in environ_only
missing_vars <- setdiff(champion_vars, names(environ_only))
if (length(missing_vars) > 0) {
  stop("The following champion variables were not found in the env table: ", paste(missing_vars, collapse = ", "))
}

# 6. Pairwise distance outputs
# First, scale the raw champion environmental variables
environ_scaled <- as.data.frame(scale(environ_only[, champion_vars, drop=FALSE]))
distmatrices <- map(environ_scaled, ~ as.matrix(dist(.x, method = "euclidean")))

# Geographic distance in Km (terra::distance on lon/lat returns meters)
geog_km <- as.matrix(terra::distance(geog)) / 1000
distmatrices <- c(distmatrices, list(geog = geog_km))

distmatrices <- map(
  distmatrices,
  ~ `dimnames<-`(.x, list(eco_var_table$Pop, eco_var_table$Pop))
)

dist_dfs_unique <- map(distmatrices, pairwise_df)

all_distances <- imap_dfr(
  dist_dfs_unique,
  ~ mutate(.x, distance_type = .y)
)

selected_vars_dist <- all_distances %>%
  pivot_wider(
    names_from  = distance_type,
    values_from = distance
  ) %>%
  arrange(pop1, pop2)

selected_vars_dist <- selected_vars_dist[, c('pop1', 'pop2', 'geog', champion_vars)]
write.csv(selected_vars_dist, args$output, row.names=FALSE, quote=FALSE)
cat("Successfully wrote scaled pairwise distances to", args$output, "\n")

# 7. Write out BayeScEnv dataset
if (!is.null(args$out_dir_txt)) {
  if (!dir.exists(args$out_dir_txt)) {
    dir.create(args$out_dir_txt, recursive = TRUE, showWarnings = FALSE)
  }
  
  final_glmm_data <- eco_var_table[, c("Pop", "lon", "lat", champion_vars)]
  

  final_glmm_scaled <- final_glmm_data %>%
    mutate(
      across(
        .cols = all_of(champion_vars),
        .fns = ~ {
          ref_val <- if (grepl("wc2_elev", cur_column())) { # Altitude variables
            0                                               # Reference = sea level
          } else {
            mean(.x, na.rm = TRUE)                          # Reference = population mean
          }
          # Compute signed contrast and standardize
          as.numeric((.x - ref_val) / sd(.x, na.rm = TRUE))
        }
      )
    )

  for (v in champion_vars) {
    write(
      x = final_glmm_scaled[[v]],
      file = file.path(args$out_dir_txt, paste0(v, "_champion.txt")),
      ncolumns = length(final_glmm_scaled[[v]])
    )
  }
  cat("Wrote txt files for BayeScEnv into directory:", args$out_dir_txt, "\n")
}
