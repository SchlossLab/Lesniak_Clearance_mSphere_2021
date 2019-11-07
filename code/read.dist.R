###################
#
# read.dist.R
#
# reads triangle matrices into R
# 
# cloned from Marc Sze
# https://github.com/SchlossLab/Sze_PCRSeqEffects_mSphere_2019/blob/master/code/mock_beta.R
#    
#
###################


library(tidyverse)

read_dist <- function(dist_file_name){
  # read in all the data from the lower triangle (exlcude the first which is the matrix dim)
  dist_data <- scan(dist_file_name, what="character", quiet=TRUE)[-1]
  # subset the sample names (which all contain "-")
  samples <- str_subset(dist_data, "-")
  n_samples <- length(samples)
  distance_strings <- str_subset(dist_data, "\\.")

  distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
  colnames(distance_matrix) <- samples
  as.tibble(cbind(rows=samples, distance_matrix)) %>%
    gather(columns, distances, -rows) %>%
    filter(rows < columns) %>%
    arrange(columns, rows) %>%
    mutate(distances = as.numeric(distance_strings))

}