#!/usr/bin/R
## This script is meant to be run directly from the Terminal using the Rscript command. Otherwise, working directory must be set.

library(rlang)
library(devtools)
library(foreach)
library(doParallel)
library(dplyr)
library(hierfstat)
library(adegenet)
library(vcfR)

# Load Funcitons:
# Function to randomly sample inds_per_sample-1 individuals per npops populations from an x hierfstat-dataframe object
random_sample<- function(x, inds_per_sample, npops){
  pop_lim = inds_per_sample
  data<-data.frame()
  current_pop = 1
  for (pop in (1:npops)){
    pop_size = 0
    data_tmp <- data.frame()
    for (indv in (1:nrow(x))){
      if (x$pop[indv] == current_pop){
        data_tmp <- rbind(x[indv,], data_tmp)
        pop_size = pop_size+1
      }
    }
    if (pop_size >= pop_lim){
      data <- rbind(data_tmp[sample(1:nrow(data_tmp), (pop_lim-1)),], data)
    } 
    else {
      data <- rbind(data_tmp[(1:nrow(data_tmp)),], data)
    }
    current_pop = current_pop + 1
  }
  return(data)
}

# Main routine:
main<-function(){
  part_result <- matrix(ncol = 4, nrow = length(2:sampl_size_max))
  row_res = 1
  for (k in (3:(sampl_size_max+1))){
    data<-random_sample(hPaxil, k, npopus)
    data<-dplyr::arrange(data, pop)
    pwfboots<-hierfstat::boot.ppfst(data, nboot=100, quant = c(0.025,0.975), diploid = TRUE)
    part_result[row_res, 1] <- k-1
    part_result[row_res, 2] <- sum(pwfboots$ll, na.rm = T)/(((npopus*npopus)-npopus)/2)
    part_result[row_res, 3] <- sum(pwfboots$ul, na.rm = T)/(((npopus*npopus)-npopus)/2)
    fstmedi<- matrix(nrow = nrow(pwfboots$ul), ncol = ncol(pwfboots$ul))
    for (i in seq(2,ncol(pwfboots$ll))){
      for (j in seq(1,nrow(pwfboots$ll))){
        fstmedi[j,i]<-(pwfboots$ll[j,i]+pwfboots$ul[j,i])/2
      }
    }
    part_result[row_res, 4] <- sum(fstmedi, na.rm = T)/(((npopus*npopus)-npopus)/2)
    row_res = row_res + 1
  }
  return(part_result)
}

# Working dir from the root of the project
# setwd('3-PopGenStruct/Subsampling_Validation/')
print(getwd())

# Read data produced with the script `script_create_hierstatdataframe.R``:
# load('Paxil_M095_hierfstat_Pops')

## load VCF:
Paxil = read.vcfR("../../2-OutlierDetection/Paxil_M095_PutatNeutral.recode.vcf.gz")
## Convert VCF to genind object:
GI_Paxil = vcfR2genind(Paxil)
## Define and assign individuals to populations:
pop_def <- read.csv('../../Paxil_Popmap.csv')
rownames(pop_def) <- pop_def$Indv
pop_def <- pop_def[rownames(GI_Paxil@tab),]
GI_Paxil@pop <- as.factor(pop_def[,2])
## Convert genind to hierfstat-dataframe object:
hPaxil <- genind2hierfstat(dat=GI_Paxil, pop = as.numeric(GI_Paxil$pop))
rownames(hPaxil) <- rownames(GI_Paxil@tab)
pop_count <- hPaxil%>%count(pop)
sampl_size_max <- max(pop_count[,2])
npopus <- length(pop_count[,2])

# Setup parallel backend to use many processors
cores = 50
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

# Start main routine:
total <- list()
total <- foreach(iteration=1:50) %dopar% {
  main()
}
# Stop cluster
stopCluster(cl)
save(total, file =  "Paxil_M095_total_runs_sample_effect")
load("Paxil_M095_total_runs_sample_effect")

# Compile the results in one table:
results <- matrix(nrow = nrow(total[[1]]), ncol = 4)
colnames(results) <- c("n", "low_limit", "high_limit", "mean")

for (j in (1:4)){
  for (i in  (1:nrow(total[[1]]))){
    vec <- c()
    for (k in (1:length(total))){
      vec[k] <- total[[k]][i,j]
    }
    results[i,j] <- mean(vec)
  }
}

results_2 <- matrix(ncol = 2, nrow = nrow(total[[1]]))
colnames(results_2) <- c('n', '95% CI')

for (i in (1:nrow(total[[1]]))){
  results_2[i,1] <- results[i,1]
  results_2[i,2] <- paste(round(results[i,4],digits = 4), '+/-', 
                        round(results[i,3]-results[i,4], digits = 4), sep = ' ')
}

write.csv(results_2, file = 'Paxil_M095_random_sample_effect.csv', row.names = F, quote = F)
