library(dartR.popgen)
library(ggplot2)
library(patchwork)

paxil_info <- read.csv('../../Paxil_PopGroupMap.csv')

hier_cols <- setdiff(colnames(paxil_info), "Indv")

paxil_pop <- unique(paxil_info[,c('Pop')])
paxil_group <- unique(paxil_info[,c('Group')])

paxil_pop_groups <- c(paxil_pop, paxil_group, "All")

plot_stairway_res <- function(file_path, prefix){
  stairway_results <- read.table(file, header=1)
  p <- ggplot(stairway_results, aes(x = year, y = Ne_median)) +
    geom_line(color = "black", linewidth=2) +
    geom_ribbon(aes(ymin = Ne_2.5., ymax = Ne_97.5.), alpha = 0.3, fill = "gray2") +
    labs(x = "Time (Years Ago)", y = "Effective Population Size (Ne)", title=prefix) +
    theme_gray(base_size=16, base_line_size=0.5, base_rect_size=0.5)
  
  return(p)
}

sum_files <- Sys.glob('results_noMono_inclPaxil1001_projCONSERVAT/*/*final.summary')
sum_files

res_plots <- list()
for (file in sum_files){
  pop_prefix <- paste(stringr::str_split(stringr::str_split_i(file,'/',2),'_')[[1]][c(6,7)], collapse='_')
  res_plots[[pop_prefix]] <- plot_stairway_res(file, pop_prefix)
}

names(res_plots)

main_stairway_plot <- (res_plots[["proj118_NA"]] | res_plots[["Cluster_A"]] | res_plots[["Cluster_B"]]) / (res_plots[["Cluster_C"]] | res_plots[["Cluster_D"]] |plot_spacer())
ggsave(filename=paste0("results_noMono_inclPaxil1001_projCONSERVAT/mainStairway_paxil_noMono_inclPaxil1001.svg"), plot=main_stairway_plot, width=15, height=10)
suppl_stairway_plot_1 <- (res_plots[["P_01"]] | res_plots[["P_02"]] | res_plots[["P_03"]]) / (res_plots[["P_04"]] | res_plots[["P_05"]] | res_plots[["P_06"]])
ggsave(filename=paste0("results_noMono_inclPaxil1001_projCONSERVAT/Suppl1Stairway_paxil_noMono_inclPaxil1001.svg"), plot=suppl_stairway_plot_1, width=15, height=10)
suppl_stairway_plot_2 <- (res_plots[["P_07"]] | res_plots[["P_08"]] | res_plots[["P_09"]]) / (res_plots[["P_10"]] | res_plots[["P_11"]] | res_plots[["P_12"]])
ggsave(filename=paste0("results_noMono_inclPaxil1001_projCONSERVAT/Suppl2Stairway_paxil_noMono_inclPaxil1001.svg"), plot=suppl_stairway_plot_2, width=15, height=10)
suppl_stairway_plot_3 <- (res_plots[["P_13"]] | plot_spacer() | plot_spacer()) / (plot_spacer() | plot_spacer() | plot_spacer())
ggsave(filename=paste0("results_noMono_inclPaxil1001_projCONSERVAT/Suppl3Stairway_paxil_noMono_inclPaxil1001.svg"), plot=suppl_stairway_plot_3, width=15, height=10)


