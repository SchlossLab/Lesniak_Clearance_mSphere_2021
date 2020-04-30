######################################################################
# Date: 04-29-2020
# Title: ML pipeline adapted from Topcuoglu mBio 2020
######################################################################

######################################################################
# Description:

# This script will read in data 
#     - Subsampled OTU data - 
#     - Metadata - 


# It will run the following machine learning pipelines:
#     - L2 Logistic Regression
######################################################################

######################################################################
# Dependencies and Outputs:

# Be in the project directory.

# The outputs are:
#   (1) AUC values for cross-validation and testing for each data-split
#   (2) meanAUC values for each hyper-parameter tested during each split.
######################################################################


################### IMPORT LIBRARIES and FUNCTIONS ###################
# The dependinces for this script are consolidated in the first part
deps = c("dplyr", "tictoc", "caret" ,"rpart", "xgboost", "randomForest", "kernlab","LiblineaR", "pROC", "tidyverse");
for (dep in deps){
	if (dep %in% installed.packages()[,"Package"] == FALSE){
		install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org", dependencies=TRUE);
	}
	library(dep, verbose=FALSE, character.only=TRUE)
}
# Load in needed functions and libraries
source('code/learning/model_selection.R')
source('code/learning/model_pipeline.R') # has pipeline function defined here (called within generate.AUCs.R)
source('code/learning/generateAUCs.R') # has get_results function defined here
source('code/learning/permutation_importance.R')
source('code/learning/compute_correlation_matrix.R')
######################################################################

######################## DATA PREPARATION #############################
# Features: 16S rRNA gene sequences(OTUs) in the stool
# Labels: - C. difficile clearance


# Read in metadata
meta <- read_tsv('data/process/abx_cdiff_metadata_clean.txt',
	col_types = 'cd--dcdlll-cdd---cddcc')  
# Read in OTU table and remove label and numOtus columns
shared <- read_tsv('data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared',
	col_types = cols(.default = 'd', Group = 'c', label = '-', numOtus = '-'))
# Filter metadata and select only sample names and clearance columns
# Merge metadata and OTU table.
# Then remove the group column
data <- meta %>% 
	filter(clearance %in% c('Cleared', 'Colonized'),
		day == 0, cdiff == T) %>% 
	select(dx = clearance, Group = group) %>% 
	inner_join(shared, by=c("Group"))  %>% 
	select(-Group) %>% 
	drop_na()

# create correlation matrix of data
#calc_corr_matrix(select(data, -dx))

# We want the diagnosis column to be a factor
data$dx <- factor(data$dx)
# We want the first sample to be a cancer so we shuffle the dataset with a specific seed to get cancer as the first sample
set.seed(0)
data <- data[sample(1:nrow(data)), ]
###################################################################

######################## RUN PIPELINE #############################
# Choose which classification methods we want to run on command line
#                "L2_Logistic_Regression",

# We will run main.R from command line with arguments
#  - These arguments will be saved into variable "input"
#  - First argument is the seed number which is the array index

input <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(input[1])
model <- 'L2_Logistic_Regression'

# Then arguments 1 and 2 will be placed respectively into the functions:
#   1. set.seed() : creates reproducibility and variability
#   2. get_results(): self-defined function that
#                     - runs the modeling pipeline
#                     - saves performance and hyper-parameters and imp features
set.seed(seed)
# Start walltime for running model
tic("model")
# Run the model
get_results(data, model, seed)
# Stop walltime for running model
secs <- toc()
# Save elapsed time
walltime <- secs$toc-secs$tic
# Save wall-time
write.csv(walltime, file=paste0("data/temp/walltime_", model, "_", seed, ".csv"), row.names=F)
###################################################################
