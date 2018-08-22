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
	rename(C_difficile = CFU)
		
# remove otus that are present in less than 10 samples
ccm_df <- select(ccm_df, day, C_difficile, which(apply(ccm_df > 1, 2, sum) > 10 ))  

normalize <- function(x, ...) {
    (x - mean(x, ...))/sd(x, ...)
}

ccm_df <- ccm_df %>% 
	# replace all 0s with random value between 0 and 1
	gather(taxa, abundance, C_difficile, contains('Otu00')) %>% 
	mutate(abundance = ifelse(abundance == 0, 
		sample(100, sum(abundance == 0), replace = T)/100, abundance)) %>% 
	spread(taxa, abundance) %>% 
	# add missing days
	full_join(data.frame(
		unique(select(ccm_df, unique_id, treatment))[
			rep(1:nrow(unique(select(ccm_df, unique_id, treatment))), each = 11), ],
		day = rep(seq(0,10), length(unique(ccm_df$unique_id))),
		stringsAsFactors = F), by = c('unique_id', 'day', 'treatment')) %>% 
	# impute values for missing days and normalize by mouse
	arrange(unique_id, day) %>% 
	gather(bacteria, raw_abundance, C_difficile, contains('Otu00')) %>% 
	group_by(unique_id, bacteria) %>% 
	nest() %>% 
	mutate(abundance = map(data, ~normalize(na.interp(.$raw_abundance)))) %>% 
	unnest(data, abundance) %>% 
	select(treatment, unique_id, day, bacteria, abundance) 

ccm_df <- ccm_df %>% 
	arrange(unique_id, bacteria, day) %>% 
	group_by(unique_id, bacteria) %>% 
	mutate(first = abundance - lag(abundance),
		second = abundance - lag(abundance, 2)) %>% 
	ungroup %>% 
	rename(raw = abundance) %>% 
	gather(differenced, abundance, raw, first, second) %>% 
	mutate(variable = paste0(bacteria, '_', differenced))

write.table(ccm_df, paste0(save_dir, 'ccm_otu_data.txt'), 
	quote = F, row.names = F)
