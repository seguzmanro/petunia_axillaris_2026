############################### 3D PCA ANALYSIS ################################

setwd("3-PopGenStruct/PCA/")

library("readr")
library("scatterplot3d")
library("factoextra")

# Read distance matrix obtained using TASSEL
Paxil_matrix_csv <- read.csv("Paxil_M095_PutatNeutral.DistMtrx.csv", row.names = 1)

Paxil_mtx = as.matrix(Paxil_matrix_csv)
Paxil_pr = prcomp(Paxil_mtx)
Paxil_pr
Paxil_pr$sdev
format(Paxil_pr$sdev * 100, digits=2, scientific = FALSE)
Paxil_pr$rotation
Paxil_pr_Rotation = Paxil_pr$rotation

eig<-get_eig(Paxil_pr)
eig$variance.percent[1]
eig$variance.percent[2]
eig$variance.percent[3]

## Prepare the data that I'll use for 3D Plot (3 most important component in PCA)
Paxil_Data = Paxil_pr_Rotation[,c(1:3)]
write.csv(x = Paxil_Data, file="Paxil_M095_PutatNeutral_Matrix_3D.csv")

## Import .csv that I've created using "Paxil_Data" (above)
popmap <- read.csv('../../Paxil_Popmap.csv')
Paxil_Matrix_3D <- read.csv("Paxil_M095_PutatNeutral_Matrix_3D.csv")
View(Paxil_Matrix_3D) 
summary(Paxil_Matrix_3D)
Paxil_Matrix_3D_pops <- merge(Paxil_Matrix_3D, popmap, by.x = 'X', by.y='Indv')
Paxil_Matrix_3D_pops$Pop <- as.factor(Paxil_Matrix_3D_pops$Pop)
Paxil_3DPlot = Paxil_Matrix_3D_pops[,c(2:5)] 

################################################################################

color_map_pops <- c(
  P_04=col2hex("goldenrod"),
  P_06=col2hex("green"), 
  P_05=col2hex("mediumorchid3"),
  P_13=col2hex("blue"),
  P_11=col2hex("violetred1"),
  P_10=col2hex("darkgray"),
  P_12=col2hex("red"),
  P_08=col2hex("yellow"),
  P_07=col2hex("lightcyan2"),
  P_09=col2hex("pink"),
  P_03=col2hex("yellowgreen"),
  P_02=col2hex("sienna2"),
  P_01="#33FFFF" 
)
point_colors <- color_map_pops[as.character(Paxil_Matrix_3D_pops$Pop)]


## Plot

colors <- as.numeric(Paxil_Matrix_3D_pops$Pop)
palette <- rainbow(length(unique(Paxil_Matrix_3D_pops$Pop)))
svg('Paxil_M095_PutatNeutral_PCA.svg', height = 6, width = 6)
scatterplot3d(Paxil_Matrix_3D_pops$PC1, Paxil_Matrix_3D_pops$PC2, Paxil_Matrix_3D_pops$PC3,
              color = point_colors,
              pch = 16, type='h',
              grid = TRUE,
              box = FALSE,
              main = "3D PCA: Outlier Markers", 
              xlab = paste0("PC1 (", round(eig$variance.percent[1], 2), "%)"),
              ylab = paste0("PC2 (", round(eig$variance.percent[2], 2), "%)"),
              zlab = paste0("PC3 (", round(eig$variance.percent[3], 2), "%)"))
sorted_names <- sort(names(color_map_pops))
sorted_colors <- color_map_pops[sorted_names]

# Add alphabetically ordered legend
legend("topright",
       legend = sorted_names,
       col = sorted_colors,
       pch = 19,
       cex = 0.8)
dev.off()
 
################################################################################