library(vcfR)
library(dplyr)

setwd('2-OutlierDetection/')

vcf_noLD <- read.vcfR('../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf.gz')

pcadapt_outliers <- read.table('PCAdapt/Paxil_M095_noLD_Outliers_PCAdapt.txt')
colnames(pcadapt_outliers) <- c('CHROM', 'POS')
bayescan_outliers <- read.table('Bayescan/Paxil_M095_noLD_Bayescan_Outliers_BA_01.txt')
colnames(bayescan_outliers) <- c('CHROM', 'POS')

common <- pcadapt_outliers[apply(pcadapt_outliers, 1, paste, collapse = "|") %in%
              apply(bayescan_outliers, 1, paste, collapse = "|"), ]

print(common)

total_outliers <- rbind(pcadapt_outliers, bayescan_outliers)
write.table(total_outliers, 'Paxil_M095_noLD_TotalOutliers.txt', sep='\t', col.names=F, row.names=F, quote=F)
