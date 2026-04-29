library(vcfR)
library(dplyr)

setwd('petunia_axillaris_2026/2-OutlierDetection')

vcf_noLD <- read.vcfR('../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf.gz')

pcadapt_outliers <- read.table('PCAdapt/Paxil_M095_noLD_Outliers_PCAdapt.txt')
colnames(pcadapt_outliers) <- c('CHROM', 'POS')
bayescan_outliers <- read.table('Bayescan/Paxil_M095_noLD_Bayescan_Outliers_BA_01.txt')
colnames(bayescan_outliers) <- c('CHROM', 'POS')

bayescenv_result_files <- list.files('../5-GenotypeEnvironAssoc/Bayescenv_champion/Results/', pattern='*Outliers_BA_01.txt', full.names = TRUE)
bayescenv_outliers <- list()
for (file in bayescenv_result_files){
  env_prefix <- gsub("../5-GenotypeEnvironAssoc/Bayescenv_champion/Results//", "", gsub("_champion_results_Outliers_BA_01.txt","",file))
  bayescenv_outliers[[env_prefix]] <- read.table(file)
  colnames(bayescenv_outliers[[env_prefix]]) <- c('CHROM', 'POS')
}

bayescenv_all <- unique(do.call(rbind, bayescenv_outliers))
write.table(bayescenv_all, 'Paxil_M095_noLD_BayeScEnvOutliers.txt', sep='\t', col.names=F, row.names=F, quote=F)

common <- pcadapt_outliers[apply(pcadapt_outliers, 1, paste, collapse = "|") %in%
              apply(bayescan_outliers, 1, paste, collapse = "|"), ]

common_bayes <- Reduce(function(x, y) merge(x, y, by = c("CHROM", "POS")),
                 list(bayescan_outliers, bayescenv_all))

print(common_bayes)

total_outliers <- rbind(pcadapt_outliers, bayescan_outliers)
write.table(total_outliers, 'Paxil_M095_noLD_TotalOutliers.txt', sep='\t', col.names=F, row.names=F, quote=F)

