#!/usr/bin/env Rscript

# Function to load or install required packages
ver_load_packages <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE, quietly = TRUE))
  need <- libs[req == FALSE]
  if(length(need) > 0) {
    install.packages(need, repos = "http://cran.us.r-project.org")
    lapply(need, require, character.only = TRUE, quietly = TRUE)
  }
}

ver_load_packages(c("argparse"))

parser <- ArgumentParser(description="Plot conStruct cross-validation results")
parser$add_argument('--sp', required=TRUE, help='Path to spatial cross-validation results txt')
parser$add_argument('--nsp', required=TRUE, help='Path to non-spatial cross-validation results txt')
parser$add_argument('--out', required=TRUE, help='Output SVG file')
parser$add_argument('--k_values', required=TRUE, nargs='+', type='integer', help='K values used in the analysis')

args <- parser$parse_args()

k_values <- args$k_values

sp.results <- as.matrix(read.table(args$sp, header = TRUE, stringsAsFactors = FALSE))
nsp.results <- as.matrix(read.table(args$nsp, header = TRUE, stringsAsFactors = FALSE))

sp.CIs <- apply(sp.results, 1, function(x){mean(x) + c(-1.96, 1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results, 1, function(x){mean(x) + c(-1.96, 1.96) * sd(x)/length(x)})

svg(args$out, height = 5, width = 10)
par(mfrow=c(1,2))

# Plot 1: All values of K
plot(rowMeans(sp.results),
     pch=19, col="blue",
     ylab="predictive accuracy", xlab="values of K",
     ylim=range(sp.results, nsp.results),
     xaxt="n",
     main="Cross-validation results (All K)")
axis(side=1, at=1:length(k_values), labels=k_values)
points(rowMeans(nsp.results), col="green", pch=19)
legend("bottomright", legend=c("Spatial", "Non-Spatial"), col=c("blue", "green"), pch=19)

# Plot 2: Zoomed in (excluding K=1 and K=2, or K=2 only if K=1 not present)
n_k <- length(k_values)
has_k1 <- 1 %in% k_values
has_k2 <- 2 %in% k_values
exclude <- c()
if (has_k2) {
    exclude <- c(exclude, which(k_values == 2))
    if (has_k1) {
        exclude <- c(exclude, which(k_values == 1))
    }
}
idx <- setdiff(1:n_k, exclude)

if (length(idx) > 1) {
    plot(rowMeans(sp.results[idx, ]),
         pch=19, col="blue",
         ylab="predictive accuracy", xlab="values of K",
         ylim=range(sp.CIs[, idx], nsp.CIs[, idx]),
         xaxt="n",
         main=paste0("Cross-validation results (K > ", min(k_values[idx]) - 1, ")"))
         
    segments(x0 = 1:length(idx),
             y0 = sp.CIs[1, idx],
             x1 = 1:length(idx),
             y1 = sp.CIs[2, idx],
             col = "blue", lwd=2)
             
    points(rowMeans(nsp.results[idx, ]), col="green", pch=19)
    
    segments(x0 = 1:length(idx),
             y0 = nsp.CIs[1, idx],
             x1 = 1:length(idx),
             y1 = nsp.CIs[2, idx],
             col = "green", lwd=2)
             
    axis(side=1, at=1:length(idx), labels=k_values[idx])
} else {
    plot(1, type="n", axes=FALSE, xlab="", ylab="", main="Zoomed plot skipped (too few K values after exclusion)")
}

dev.off()
