library(tidyverse)
library(forecast)

save_dir <- paste0('data/process/')

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(group, cage, mouse, day, CFU, cdiff, abx, dose, delayed) %>% 
	unite(treatment, abx, dose, delayed) %>% 
	filter(cdiff == T, day >= 0, treatment != 'none_NA_FALSE')
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F)

# normalize each taxa for each subject
normalize <- function(x, ...) {
	(x - mean(x, ...))/sd(x, ...)
}

ccm_df <- meta_file %>% 
	mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
	inner_join(shared_file, by = c('group' = "Group")) %>% 
	select(-group, -numOtus) %>%
	rename(C_difficile = CFU) %>%  
	gather(taxa, raw_abundance, C_difficile, contains('Otu00')) %>% 
	group_by(taxa) %>% 
	mutate(normalized_abundance = normalize(raw_abundance)) 

taxa_by_treatment <- ccm_df %>% 
	group_by(treatment, taxa) %>% 
	# remove otus that are present in less than 10 samples for each treatment	
	mutate(present = sum(raw_abundance > 1) > 10 ) %>% 
	ungroup() %>% 
	filter(present == T) %>% 
	# fill in missing days
	select(unique_id, treatment, taxa, day) %>% 
	mutate(present = T) %>% 
	spread(day, present) %>% 
	gather(day, present, -unique_id, -treatment, -taxa) %>% 
	select(-present) %>% 
	mutate(day = as.numeric(day))

ccm_df <- ccm_df %>% 
	# add missing days and only use specific otus per treatment group
	right_join(taxa_by_treatment, by = c('unique_id', 'day', 'treatment', 'taxa')) %>% 
	# impute values for missing days and normalize by mouse
	arrange(unique_id, taxa, day) %>% 
	group_by(unique_id, taxa) %>% 
	nest() %>% 
	mutate(normalized_abundance = map(data, ~na.interp(.$normalized_abundance))) %>% 
	unnest(data, normalized_abundance) %>% 
	select(treatment, unique_id, day, taxa, normalized_abundance, raw_abundance)

ccm_df <- ccm_df %>% 
	arrange(unique_id, taxa, day) %>% 
	group_by(unique_id, taxa) %>% 
	mutate(first = normalized_abundance - lag(normalized_abundance),
		second = normalized_abundance - lag(normalized_abundance, 2)) %>% 
	ungroup %>% 
	rename(undifferenced = normalized_abundance) %>% 
	gather(differenced, normalized_abundance, undifferenced, first, second) %>% 
	mutate(otu_feature = paste0(taxa, '_', differenced))

write.table(ccm_df, paste0(save_dir, 'ccm_otu_data.txt'), 
	quote = F, row.names = F)
