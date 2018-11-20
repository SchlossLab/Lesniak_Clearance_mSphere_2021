
# load interaction data
interactions_df <- c()
for(i in unique(xmap_otus$driven)){ 
	tmp <- read.table(paste0(save_dir, '/interactions_w_', i, '.txt'), header = T)
	tmp <- tmp %>% filter(embed == 'multi') %>% 
	gather(interaction, strength, contains('dOTU')) %>% 
	separate(interaction, c('junk1', 'affected', 'junk2','affector')) %>% 
	select(affected, affector, strength) %>% 
	mutate(affected = paste0('OTU_', affected),
		affector = paste0('OTU_', affector)) %>% 
	inner_join(true_interactions, by = c('affected' = 'affected_otu', 'affector' = 'interaction'))  %>% 
	filter(!is.na(actual_strength), !is.na(strength))
	interactions_df <- bind_rows(interactions_df, tmp)
}

interactions_df %>% 
	ggplot() + 
	geom_point(aes(actual_strength, strength), alpha = 0.01)