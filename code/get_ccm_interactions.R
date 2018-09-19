library(rEDM)
library(tidyverse)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
print(paste0('Running set ', run_set))

otu_df <- 'data/process/ccm_otu_data.txt'
otu_df   <- read.table(otu_df, header = T, stringsAsFactors = F) %>% 
	mutate(otu_feature = gsub('_first', '', otu_feature))

seed <- 062818
treatment_subset <- unique(otu_df$treatment)[run_set]
print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset))
otu_df <- filter(otu_df, treatment %in% treatment_subset)
save_dir <- paste0('scratch/ccm_otu/', treatment_subset, '/interactions')
ifelse(!dir.exists(save_dir), 
	dir.create(save_dir), 
	print(paste0(save_dir, ' directory ready')))

# import results from individual taxa test
#	check for embedding/nonlinearity/ccm files
load_files <- c('simplex_embedding_first_differenced.txt',
	'smap_nonlinearity_first_differenced.txt',
	paste0('ccm_by_otu_', treatment_subset, '_first_differenced.txt'))
for(i in load_files){
	if(!file.exists(paste0(save_dir, '/../', i))){ 
		stop(paste('Run previous script or locate', i))
		}
	}
#  import embedding/nonlinearity
nonlinearity <- read.table(paste0(save_dir, '/../', load_files[2]), 
	header = T, stringsAsFactors = F) %>% 
		mutate(taxa = gsub('_first', '', taxa))
#  import ccm analysis
ccm <- read.table(paste0(save_dir, '/../', load_files[3]), 
	header = T, stringsAsFactors = F) 

# check for significance in all tests
#	nonlinear
#		decrease in mae (increase in prediciton with increased theta) (wilcox)
nonlinear_otus <- nonlinearity %>% 
	filter(p_linear_v_nonlinear < 0.05, 
		p_real_v_surrogate < 0.05, 
		data == 'real') %>% 
		pull(taxa)

#	ccm convergence
#		increasing (spearman), full library > min library (wilcox)
#	ccm significance
#		greater than linear (ttest), greater than null model (wilcox)
# confirm median of real data is higher
xmap_otus <- ccm %>% 
	filter(data == 'real',
		ccm_trend_p < 0.05,
		ccm_trend_rho > 0,
		median_rho_max_lobs > 0.25,
		p_min_v_max_lobs < 0.05,
		linear_corr_p < 0.05,
		ccm_null_p < 0.05) %>% 
	select(causal) %>% 
	separate(causal, sep = 'xmap', c('driven', 'driver'))

# for all interactions that are significant
# import embeddings
# check mae for best theta
# run smap with best theta
# partial derivaties of multivariate smap approximates interactions

# output file with taxa, interaction direction, interaction strength

