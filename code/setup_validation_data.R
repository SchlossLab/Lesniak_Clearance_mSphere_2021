library(tidyverse)

save_dir <- paste0('data/process/')

shared_file <- 'data/process/validation/validation_temporal_data_gamma2_seed3.txt'
shared_file <- read.table(shared_file, stringsAsFactors = F, header = T)

shared_file <- shared_file %>% 
	mutate(sample_id = 1:NROW(shared_file))

# normalize each taxa for each subject
normalize <- function(x, ...) {
	(x - mean(x, ...))/sd(x, ...)
}

ccm_df <- shared_file %>%  
	rename(unique_id = replicate, day = time) %>% 
	gather(taxa, raw_abundance, contains('OTU')) %>% 
	group_by(taxa) %>% 
	mutate(normalized_abundance = normalize(raw_abundance))  %>% 
	# remove otus that are present in less than 10 samples for each treatment	
	mutate(present = sum(raw_abundance > 1) > 10 ) %>% 
	ungroup() %>% 
	filter(present == T)
# test missing days
						#%>% 
#	# fill in missing days
#	select(unique_id, taxa, day) %>% 
#	mutate(present = T) %>% 
#	spread(day, present) %>% 
#	gather(day, present, -unique_id, -taxa) %>% 
#	select(-present) %>% 
#	mutate(day = as.numeric(day))

#ccm_df <- ccm_df %>% 
#	# add missing days and only use specific otus per treatment group
#	right_join(taxa_by_treatment, by = c('unique_id', 'day', 'taxa')) %>% 
#	# impute values for missing days and normalize by mouse
#	arrange(unique_id, taxa, day) %>% 
#	group_by(unique_id, taxa) %>% 
#	nest() %>% 
#	mutate(normalized_abundance = map(data, ~na.interp(.$normalized_abundance))) %>% 
#	unnest(data, normalized_abundance) %>% 
#	select(treatment, unique_id, day, taxa, normalized_abundance, raw_abundance)

ccm_df <- ccm_df %>% 
	arrange(unique_id, taxa, day) %>% 
	group_by(unique_id, taxa) %>% 
	mutate(first = normalized_abundance - lag(normalized_abundance)) %>% 
	ungroup %>% 
	rename(undifferenced = normalized_abundance) %>% 
	gather(differenced, normalized_abundance, undifferenced, first) %>% 
	mutate(otu_feature = paste0(taxa, '_', differenced),
		treatment = 'validation')

write.table(ccm_df, paste0(save_dir, 'ccm_validation_data.txt'), 
	quote = F, row.names = F)
