library(dplyr)
library(tidyr)
library(fuzzySim)
library(purrr)
library(terra)

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

setwd('5-GenotypeEnvironAssoc/EnvironValriables/')

eco_var_table <- read.csv(
  "Paxil_EnvironVars.csv",
  stringsAsFactors = TRUE
)

paxil_popmap <- read.csv('../../Paxil_PopGroupMap.csv')

rownames(eco_var_table) <- eco_var_table$Pop
eco_var_table <- eco_var_table[paxil_popmap[,'Pop'][!duplicated(paxil_popmap[,'Pop'])],]


geog <- vect(
  eco_var_table[, c("lon", "lat")],
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

environ_only <- eco_var_table[, grepl("mean|CM|wc|Alt", names(eco_var_table))]

distmatrices <- c(
  # Environmental Euclidean distances
  map(environ_only, ~ as.matrix(dist(.x, method = "euclidean"))),

  # Geographic great-circle distances
  list(geog = as.matrix(terra::distance(geog)))
)

distmatrices <- map(
  distmatrices,
  ~ `dimnames<-`(.x, list(eco_var_table$Pop, eco_var_table$Pop))
)

dist_dfs_unique <- map(distmatrices, pairwise_df)

all_distances <- imap_dfr(
  dist_dfs_unique,
  ~ mutate(.x, distance_type = .y)
)

dist_wide <- all_distances %>%
  pivot_wider(
    names_from  = distance_type,
    values_from = distance
  ) %>%
  arrange(pop1, pop2)

dist_scaled <- dist_wide %>%
  mutate(
    across(
      -c(pop1, pop2),
      ~ as.numeric(scale(.x))
    )
  )



soil_pairwise_dist <- dist_scaled[,grepl('mean|geog',names(dist_scaled))]

soil_selec_vars_VIF <- corSelect(soil_pairwise_dist, 
          sp.cols = length(names(soil_pairwise_dist)), 
          var.cols = 1:(length(names(soil_pairwise_dist))-1), 
          method='pearson', select='VIF', coeff = FALSE, cor.thresh=0.01)

soil_selec_vars_VIF$remaining.multicollinearity


atmos_pairwise_dist <- dist_scaled[,grepl('CM|geog',names(dist_scaled))]
atmos_selec_vars_VIF <- corSelect(atmos_pairwise_dist, 
                                              sp.cols = length(names(atmos_pairwise_dist)), 
                                              var.cols = 1:(length(names(atmos_pairwise_dist))-1), 
                                              method='pearson', select='VIF', coeff = FALSE,cor.thresh=0.01)

atmos_selec_vars_VIF$remaining.multicollinearity


selected_vars_dist <- dist_scaled[c(
    'pop1','pop2','geog',
    'wc2_elev','Alt_tri','Alt_topoWet',
    atmos_selec_vars_VIF$selected.vars, 
    soil_selec_vars_VIF$selected.vars
)]

selected_vars_raw <- eco_var_table[c(
    'Pop', 'lon','lat',
    'wc2_elev','Alt_tri','Alt_topoWet',
    atmos_selec_vars_VIF$selected.vars, 
    soil_selec_vars_VIF$selected.vars
)]

write.csv(selected_vars_dist, 'EnvironVars_Selected_Distances.csv', row.names=FALSE, quote=FALSE)

env_dist_to_mean <- selected_vars_raw %>%
  mutate(
    across(
      -c(Pop, lon, lat,),
      ~ (.x - mean(.x, na.rm = TRUE))
    )
  )

env_dist_to_mean <- env_dist_to_mean %>%
  mutate(Pop = selected_vars_raw$Pop) %>%
  relocate(Pop)

env_dist_scaled <- env_dist_to_mean %>% 
  mutate(
    across(
      -c(Pop, lon, lat),
      ~ as.numeric(abs(scale(.x)))
    )
  )

vars_to_write <- setdiff(
  names(env_dist_scaled),
  c("Pop", "lon", "lat")
)

for (v in vars_to_write) {
  write(
    x = env_dist_scaled[[v]],
    file = paste0('../Bayescenv/', v, ".txt"),
    ncolumns = length(env_dist_scaled[[v]])
  )
}
