library(pcadapt)
library(fsthet)
library(qqman)
library(qvalue)
library(vcfR)
library(dplyr)


setwd('2-OutlierDetection/PCAdapt/')

# Read VCF
paxil_M095_vcf = read.vcfR('../../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf')
paxil_M095_pcadapt = read.pcadapt('../../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf', type='vcf')

# Make PCAdapt with 60 Principal Components (PCs) (one can use as many as one wants, usually people use only around 30 or so PCs); the min.maf parameter should be set equal to the maf used when filtering the VCF:
PaxilOnly_Nclusters = pcadapt(input = paxil_M095_pcadapt, K = 60, min.maf=0.01)
# Make a screeplot to choose the number of PCs using the "elbow rule"
plot(PaxilOnly_Nclusters, option = "screeplot", K = 60)

# Population map in CSV, with header. The first column has the Individual name (as it appears in the VCF) and the second the Population it belongs to.
paxil_info <- read.csv('../../Paxil_Popmap.csv')
colnames(paxil_info) <- c('indv','pop')

# Order the PopMap according to the VCF sample names
row.names(paxil_info) <- paxil_info$indv
paxil_info <- paxil_info[colnames(paxil_M095_vcf@gt)[-c(1)],]

# Plot the first 5 PCs to check they are indeed significant for our samples:
plot(PaxilOnly_Nclusters, option = "scores", i=2, j=1, pop = paxil_info$pop) 
plot(PaxilOnly_Nclusters, option = "scores", i=2, j=3, pop = paxil_info$pop) 
plot(PaxilOnly_Nclusters, option = "scores", i=4, j=3, pop = paxil_info$pop) 
plot(PaxilOnly_Nclusters, option = "scores", i=4, j=5, pop = paxil_info$pop) 


PaxilOnly_K13 = pcadapt(input = paxil_M095_pcadapt, K = 13, min.maf=0.01)
summary(PaxilOnly_K13)

K13_mtt_plot = plot(PaxilOnly_K13 , option = "manhattan")
data_K13_mtt_plot = K13_mtt_plot$data

# QQ plot
plot(PaxilOnly_K13, option = "qqplot")

# Histogram of the distribution of p.values
hist(PaxilOnly_K13$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "orange")

alpha = 0.01
# qvalue correction
qval = qvalue(PaxilOnly_K13$pvalues)$qvalues
outliers_q_01 = which(qval < alpha)
length(outliers_q_01)
# Benjamini-Hochberg correction
PaxilOnlyj_K13_bh = p.adjust(PaxilOnly_K13$pvalues,method="BH")
outliers_bh_01 = which(PaxilOnlyj_K13_bh < alpha)
length(outliers_bh_01)
# BY correction
PaxilOnlyj_K13_by = p.adjust(PaxilOnly_K13$pvalues,method="BY")
outliers_by_01_K13 = which(PaxilOnlyj_K13_by < alpha)
length(outliers_by_01_K13)
# Bonferroni correction
PaxilOnlyj_K13_bonf = p.adjust(PaxilOnly_K13$pvalues,method="bonferroni")
outliers_bonf_01 = which(PaxilOnlyj_K13_bonf < alpha)
length(outliers_bonf_01)

PaxilOnly_K13_excNA = na.exclude(PaxilOnly_K13$pvalues)
PaxilOnly_K13_manhattan = K13_mtt_plot$data$y
PaxilOnly_K13_excNA_bonf = p.adjust(PaxilOnly_K13_excNA, method="bonferroni")
PaxilOnly_K13_dat = data.frame(cbind(PaxilOnly_K13$pass, PaxilOnly_K13_manhattan, PaxilOnly_K13_excNA, PaxilOnly_K13_excNA_bonf))
names(PaxilOnly_K13_dat) =  c("snp", "p-value_manhattan", "p-value", "adjusted_p-value")
summary(PaxilOnly_K13_dat)
# In the end, as my p-values were really, REALLY small (many smaller than 1e-323), I had to filter them 
# by the log(p-value), which is computed by (and only accessible through) the manhattan-plot function and filtered 
# the loci with 'manhattan' p-value > 2000 [ log(p-value)>2000 ], which seems reasonable in the manhattan plot
data_PaxilOnly_K13_thresh = PaxilOnly_K13_dat[PaxilOnly_K13_dat$`p-value_manhattan`>2000,]
nrow(data_PaxilOnly_K13_thresh)

# Plot the manhattan plot and overlay the selected loci in red
plot(PaxilOnly_K13_dat$'p-value_manhattan'~PaxilOnly_K13_dat$snp, pch=16,
     ylim = c(0,max(PaxilOnly_K13_dat$`p-value_manhattan`)), xlim = c(0,(length(c(PaxilOnly_K13$pvalues)))), ylab = "-log10(p-values)", xlab = "SNPs")

par(new=T)

plot(data_PaxilOnly_K13_thresh$'p-value_manhattan'~data_PaxilOnly_K13_thresh$snp, pch=16, col="red",
     ylab="", xlab="", ylim = c(0,max(PaxilOnly_K13_dat$`p-value_manhattan`)), xlim = c(0,(length(c(PaxilOnly_K13$pvalues)))))

row.names(paxil_M095_vcf@fix) <- c(1:nrow(paxil_M095_vcf@fix))
write.table(paxil_M095_vcf@fix[row.names(data_PaxilOnly_K13_thresh),c(1,2)], file = 'Paxil_M095_noLD_Outliers_PCAdapt.txt', sep='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

