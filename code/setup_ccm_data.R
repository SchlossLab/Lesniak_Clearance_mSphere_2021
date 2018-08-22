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

ccm_df <- meta_file %>% 
	mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
	inner_join(shared_file, by = c('group' = "Group")) %>% 
	select(-group, -numOtus) %>% 
	rename(C_difficile = CFU) %>% 
	gather(taxa, raw_abundance, C_difficile, contains('Otu00'))

		

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

# normalize each taxa for each subject
normalize <- function(x, ...) {
    (x - mean(x, ...))/sd(x, ...)
}
ccm_df <- ccm_df %>% 
	# add missing days and only use specific otus per treatment group
	right_join(taxa_by_treatment, by = c('unique_id', 'day', 'treatment', 'taxa')) %>% 
	# impute values for missing days and normalize by mouse
	arrange(unique_id, taxa, day) %>% 
	group_by(unique_id, taxa) %>% 
	nest() %>% 
	mutate(abundance = map(data, ~normalize(na.interp(.$raw_abundance)))) %>% 
	unnest(data, abundance) %>% 
# dont think i need this since im not using mccm and now using z-score normalized data
#	# replace all 0s with random value between 0 and 1
#	mutate(abundance = ifelse(abundance == 0, 
#		sample(100, sum(abundance == 0), replace = T, na.rm = T)/100, abundance)) %>% 
	select(treatment, unique_id, day, taxa, abundance) 

ccm_df <- ccm_df %>% 
	arrange(unique_id, taxa, day) %>% 
	group_by(unique_id, taxa) %>% 
	mutate(first = abundance - lag(abundance),
		second = abundance - lag(abundance, 2)) %>% 
	ungroup %>% 
	rename(raw = abundance) %>% 
	gather(differenced, abundance, raw, first, second) %>% 
	mutate(variable = paste0(taxa, '_', differenced))

write.table(ccm_df, paste0(save_dir, 'ccm_otu_data.txt'), 
	quote = F, row.names = F,)
