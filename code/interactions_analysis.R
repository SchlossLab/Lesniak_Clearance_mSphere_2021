library(tidyverse)


input_dir <- 'data/process/ccm/interactions/'
file_list <- list.files(path = input_dir)

all_interactions <- map_dfr(file_list, function(i){
	read.table(paste0(input_dir, i), header = T, stringsAsFactors = F) %>%
		mutate(grouping = gsub('.txt', '', i)) %>%
		mutate(grouping = gsub('interactions_', '', grouping)) # %>%
#		separate(file, c('gamma', 'seed')) %>%
#		mutate(gamma = gsub('gamma', '', gamma)) %>%
#		mutate(seed = gsub('seed', '', seed))
		}) 

get_sens_spec <- function(input_data){
	data.frame(
		sample = unique(input_data$grouping),
		Actual_pos =  sum(input_data$actual_ixn_TF == T),
		Actual_neg =  sum(input_data$actual_ixn_TF == F),
		TP = sum(input_data$actual_ixn_TF == T & input_data$predicted_ixn_TF == T), 
		TN = sum(input_data$actual_ixn_TF == F & input_data$predicted_ixn_TF == F), 
		FP = sum(input_data$actual_ixn_TF == F & input_data$predicted_ixn_TF == T), 
		FN = sum(input_data$actual_ixn_TF == T & input_data$predicted_ixn_TF == F), 
		# Sensitivity - probability test will id interaction among those interacting
		# Sensitivity: TP / ( TP + FN ) × 100
		# Sensitivity: TP / ( TrueIteraction ) × 100
		# Specificity - fraction not intetacting identified not to interact
		# Specificity: TN / ( TN + FP ) × 100
		# Specificity: TN / ( True_NON_interaction ) × 100
		sensitivity = sum(input_data$actual_ixn_TF == T & input_data$predicted_ixn_TF == T) / 
	  		sum(input_data$actual_ixn_TF == T),
		specificity = sum(input_data$actual_ixn_TF == F & input_data$predicted_ixn_TF == F) / 
			sum(input_data$actual_ixn_TF == F),
		stringsAsFactors = F
	)
}

bind_rows(all_interactions %>%
		list %>%
		map_dfr(., get_sens_spec) %>%
		mutate(sample = 'all') %>%
		unique,
		bind_rows(all_interactions %>%
			split(., .$grouping) %>%
			map_dfr(., get_sens_spec), 
			all_interactions %>%
			separate(grouping, c('gamma', 'seed')) %>%
			mutate(grouping = gamma) %>%
			split(., .$grouping) %>%
			map_dfr(., get_sens_spec)
			)
		)

## load interaction data
#interactions_df <- c()
#for(i in unique(xmap_otus$driven)){ 
#	tmp <- read.table(paste0(save_dir, '/interactions_w_', i, '.txt'), header = T)
#	tmp <- tmp %>% filter(embed == 'multi') %>% 
#	gather(interaction, strength, contains('dOTU')) %>% 
#	separate(interaction, c('junk1', 'affected', 'junk2','affector')) %>% 
#	select(affected, affector, strength) %>% 
#	mutate(affected = paste0('OTU_', affected),
#		affector = paste0('OTU_', affector)) %>% 
#	inner_join(true_interactions, by = c('affected' = 'affected_otu', 'affector' = 'interaction'))  %>% 
#	filter(!is.na(actual_strength), !is.na(strength))
#	interactions_df <- bind_rows(interactions_df, tmp)
#}
#
#interactions_df %>% 
#	ggplot() + 
#	geom_point(aes(actual_strength, strength), alpha = 0.01)