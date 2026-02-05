source('~/install_load.R')

list.of.packages <- c('vcfR',
                      'hierfstat',
                      'poppr',
                      'dotenv',
                      'argparse')

install_load(list.of.packages)

if (sys.nframe() == 0) {
  
  parser <- ArgumentParser(description='Hierfstat R script for Pop Genomics - 17-10-2024')
  parser$add_argument('--vcf', type='character', help='VCF ENV unique substring')
  parser$add_argument('--popmap', type='character', help='Population map CSV. The first column must have sample names and the second column the population they belong to')
  parser$add_argument('--out_dir', type='character', help='Name of output folder')
  parser$add_argument('--out_prefix', type='character', help='Prefix of output files')
  
  args <- parser$parse_args()
  if  (!endsWith(args$out_dir, '/')){
    args$out_dir <- paste0(args$out_dir,'/')
  }
  
  print(args)
  
  loaded_vcf <- read.vcfR(args$vcf) ## load VCF
  loaded_genind <- vcfR2genind(loaded_vcf)   ## convert VCF to genind object
  pop_table <- read.csv(args$popmap)
  if (ncol(pop_table) == 3){
    colnames(pop_table) <- c('indv','pop','group')
  } else {
    colnames(pop_table) <- c('indv','pop')
  }
  
  dir.create(args$out_dir, recursive = TRUE)
  setwd(args$out_dir)
  
  table(pop_table$pop)
  pop(loaded_genind) = pop_table$pop
  
  Stats01 = basic.stats(loaded_genind)
  
  Stats01$pop.freq
  levels(loaded_genind@pop)
  # Confidence Interval for Fis:
  Fis_ci <-boot.ppfis(dat=loaded_genind, nboot=100, quant = c(0.025,0.975), diploid = TRUE, dig = 4)

  #Hs Remove Hs with 'NaN'
  Hs1 <- matrix(nrow = nrow(Stats01$Fis), ncol = ncol(Stats01$Fis))
  for (i in seq(1,ncol(Stats01$Fis))){
    for (j in seq(1,nrow(Stats01$Fis))){
      if (Stats01$Fis[j,i] != 'NaN'){
        Hs1[j,i] <- Stats01$Hs[j,i]
      } 
    }
  }
  
  # Ho Remove Ho with 'NaN'
  Ho1 <- matrix(nrow = nrow(Stats01$Fis), ncol = ncol(Stats01$Fis))
  for (i in seq(1,ncol(Stats01$Fis))){
    for (j in seq(1,nrow(Stats01$Fis))){
      if (Stats01$Fis[j,i] != 'NaN'){
        Ho1[j,i] <- Stats01$Ho[j,i]  
      } 
    }
  }
  
  # Compiling all results in one table
  Result_matrix <- matrix(nrow = length(levels(loaded_genind@pop)), ncol = 5) # nrow = number of populations 
  pop_names <- levels(loaded_genind@pop)  
  colnames(Result_matrix) <- c('Hs','Ho', 'Fis_(Man_calc)', 'Low-Lim', 'High-Lim')
  rownames(Result_matrix) <- pop_names
  for (i in 1:ncol(Stats01$Fis)){
    Result_matrix[i,2]<-mean(Ho1[,i], na.rm = TRUE);
    Result_matrix[i,1]<-mean(Hs1[,i], na.rm = TRUE)
  }
  
  for (i in seq(1,nrow(Result_matrix))){
    Result_matrix[i,3]<-1-(Result_matrix[i,2]/Result_matrix[i,1])  
    Result_matrix[i,4]<-Fis_ci$fis.ci[i,1]
    Result_matrix[i,5]<-Fis_ci$fis.ci[i,2]
  }
  
  write.csv(Stats01$overall, paste0(args$out_prefix,'_overall.csv'))
  write.csv(Result_matrix, file = paste0(args$out_prefix,'_PopStats.csv'))
  
  hierf_data <- genind2hierfstat(dat = loaded_genind, pop = as.numeric(loaded_genind$pop))
  pwfboots <- boot.ppfst(hierf_data, nboot=100, quant = c(0.025,0.975), diploid = TRUE)
  colnames(pwfboots$ll) <- levels(loaded_genind$pop)
  row.names(pwfboots$ll) <- levels(loaded_genind$pop)
  colnames(pwfboots$ul) <- levels(loaded_genind$pop)
  row.names(pwfboots$ul) <- levels(loaded_genind$pop)
  write.csv(pwfboots$ll, file = paste0(args$out_prefix,'_fst_signif_ll.csv'))
  write.csv(pwfboots$ul, file = paste0(args$out_prefix,'_fst_signif_ul.csv'))
  
  fstmedi<- matrix(nrow = nrow(pwfboots$ul), ncol = ncol(pwfboots$ul))
  colnames(fstmedi) <- levels(loaded_genind$pop)
  rownames(fstmedi) <- levels(loaded_genind$pop)
  for (i in seq(2,ncol(pwfboots$ll))){
    for (j in seq(1,nrow(pwfboots$ll))){
      fstmedi[j,i]<-(pwfboots$ll[j,i]+pwfboots$ul[j,i])/2
    }
  }
  for (i in 1:(nrow(pwfboots$ll)-1)){
    for (j in 2:ncol(pwfboots$ll)){
      if ((isTRUE(pwfboots$ll[i,j]>0) && (isTRUE(pwfboots$ul[i,j]>0)))||
          ((isTRUE(pwfboots$ll[i,j]<0)&&(isTRUE(pwfboots$ul[i,j]<0))))){
        fstmedi[j,i]<-"*"
      }
    }
  }
  
  write.csv(fstmedi, file = paste0(args$out_prefix,'_fst_signif_medi.csv'))
  # AMOVA
  if (ncol(pop_table) == 3){
    strata(loaded_genind) <- pop_table
  
    amova_res<- poppr.amova(loaded_genind, ~group/pop, nperm = 99999, method = 'pegas', within = FALSE)
    
    
    save(amova_res, file = paste0(args$out_prefix,'_amova_res'))
    write.csv(amova_res$varcomp, file=paste0(args$out_prefix, '_amova_var_components.csv'))
    
  }
}