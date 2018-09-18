library(rEDM)
library(tidyverse)
library(gtools)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
save_dir <- paste0('scratch/ccm_otu/')
print(paste0('Running set ', run_set))

# import results from individual taxa test
#  import embedding/nonlinearity
#  import ccm analysis

# check for significance in all tests
#	nonlinear
#		decrease in mae (increase in prediciton with increased theta) (wilcox)
#	ccm convergence
#		increasing (spearman), full library > min library (wilcox)
#	ccm significance
#		greater than linear (ttest), greater than null model (wilcox)
# confirm median of real data is higher

# for all interactions that are significant
# import embeddings
# check mae for best theta
# run smap with best theta
# partial derivaties of multivariate smap approximates interactions

# output file with taxa, interaction direction, interaction strength

