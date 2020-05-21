######################################################################
# Date: 05-04-2020
# Author: Nick Lesniak
# Title: Setup data to be run in machine learning pipeline
######################################################################

######################################################################
# Description:

# This script will read in data 
#     - Feature data (shared file w/OTUs)
#     - Metadata 


# It will run the following:
#     - code/R/compute_correlation_matrix.R
######################################################################

######################################################################
# Dependencies and Outputs:

# Be in the project directory.

# The outputs are:
#   (1) data/process/LEVEL_input_data.csv - CSV with data for machine learning -
#			first column is outcome of interest, 
#			remaining columns are features, one per column.
#   (2) data/process/sig_flat_corr_matrix_LEVEL.csv - CSV with correlated features
######################################################################

meta_file <- 'data/process/abx_cdiff_metadata_clean.txt'
feature_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'

level <- 'otu'

################### IMPORT LIBRARIES and FUNCTIONS ###################
# The dependinces for this script are consolidated in the first part
deps = c("tidyverse", "caret", "Hmisc");
for (dep in deps){
	if (dep %in% installed.packages()[,"Package"] == FALSE){
		install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org", dependencies=TRUE);
	}
	library(dep, verbose=FALSE, character.only=TRUE)
}
# Load in needed functions and libraries
source('code/R/compute_correlation_matrix.R')
source('code/sum_otu_by_taxa.R')
######################################################################


######################## DATA PREPARATION #############################


# ----------------------- Read in data --------------------------------
# Read in metadata
meta_df <- read_tsv(meta_file,
	col_types = 'cdc-dcdlll-cdd---cddcc')

# Read in OTU table and remove label and numOtus columns
features_df <- read_tsv(feature_file,
	col_types = cols(.default = 'd', Group = 'c', label = '-', numOtus = '-'))

# Read in taxonomy data
tax_df <- read_tsv(tax_file, col_types = cols(.default = 'c')) 

# ---------------------------------------------------------------------


# ----------------- Select samples and features -----------------------
# Filter metadata and select only sample names and outcome columns
 model_df <- meta_df %>% 
	filter(abx %in% c('Clindamycin', 'Cefoperazone', 'Streptomycin'),
		clearance %in% c('Cleared', 'Colonized'),
		day == 0, cdiff == T) %>% 
	mutate(cage = paste(grep('^.{3}', abx, value = TRUE), dose_level, 'cage', cage, sep = '_')) %>% 
	select(Group = group, cage, abx, clearance) %>% 
# Create features for potentially confounding variables
	mutate(abx_used = 1, cage_present = 1) %>% 
	spread(abx, abx_used, fill = 0) %>% 
	spread(cage, cage_present, fill = 0)
	
# Merge metadata and feature data.
# Then remove the sample name column
otu_data <- model_df %>% 
	inner_join(features_df, by = "Group")
	
#genus_data <- sum_otu_by_taxa(tax_df, features_df, 'Genus') %>% 
#	spread(taxa, abundance) %>% 
#	inner_join(model_df, by = 'Group')
#
#family_data <- sum_otu_by_taxa(tax_df, features_df, 'Family') %>% 
#	spread(taxa, abundance) %>% 
#	inner_join(model_df, by = 'Group')
#
#order_data <- sum_otu_by_taxa(tax_df, features_df, 'Order') %>% 
#	spread(taxa, abundance) %>% 
#	inner_join(model_df, by = 'Group')

# save names of row samples
otu_data %>% 
	select(Group) %>% 
	left_join(meta_df, by = c('Group' = 'group')) %>% 
	select(Group, cage, clearance) %>% 
	write_csv('data/process/', level, '_sample_names.txt') %>% 
	drop_na()
# ---------------------------------------------------------------------


save_data <- function(data, level){
	data <- data %>% 
		select(-Group)
	# ---------------------- Process model data ---------------------------
	# Remove features with near zero variance and scale remaining from 0 to 1
	preProcValues <- preProcess(data, method = c("nzv", "range"))
	dataTransformed <- predict(preProcValues, data)
	# Save data to be used in machine learning pipeline
	write_csv(dataTransformed, paste0('data/process/', level, '_input_data.csv'))
	# ---------------------------------------------------------------------


	# ------------------- Create correlation matrix -----------------------
	# Create correlation matrix of machine learning data
	#   filters correlation >= cor_value and p values < p_value
	#   default values are cor_value = 1, and p_value = 0.1
	compute_correlation_matrix(input_file = paste0('data/process/', level, '_input_data.csv'), 
		outcome = 'clearance', level = level,
		cor_value = 0.8, p_value = 0.05)
	# ---------------------------------------------------------------------
}

save_data(otu_data, level)
#save_data(genus_data, 'genus')
#save_data(family_data, 'family')
#save_data(order_data, 'order')
