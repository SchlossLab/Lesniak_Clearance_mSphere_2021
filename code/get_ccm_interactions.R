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

meta_file <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_file <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(group, cage, mouse, day, CFU, cdiff, abx, dose, delayed) %>% 
	unite(treatment, abx, dose, delayed) %>% 
	filter(cdiff == T, day >= 0, treatment == treatment_subset)
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	inner_join(meta_file, by = c("Group" = 'group')) %>%
	rename(C_difficile = CFU)

mouse_list <- unique(otu_df$unique_id)

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

for(i in unique(xmap_otus$driven)){
	otus <- xmap_otus %>% 
		filter(driven == i) %>% 
		pull(driver)
	print(paste(i, 'is driven by', paste(otus, collapse = ', ')))
	dynamics_plot <- shared_file %>% 
			select(cage, mouse, day, one_of(i, otus)) %>% 
			gather(bacteria, counts, one_of(i, otus)) %>% 
				ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), group = interaction(cage, mouse))) + 
					geom_line(alpha = 0.4) + 
					geom_point() + 
					facet_grid(bacteria~., scales = 'free_y') +
					theme_bw() + 
					labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
						subtitle = 'Temporal Dynamics - Colored by mouse', 
						title = print(paste(i, 'is driven by', paste(otus, collapse = ', ')))) + 
					scale_x_continuous(breaks=seq(0,10, 1)) + 
					theme_bw(base_size = 8) + 
					theme(legend.position = 'none', panel.grid.minor = element_blank())
	
	composite_ts <- otu_df %>% 
		filter(taxa %in% c(i, otus), differenced == 'first') %>% 
		select(taxa, day, unique_id, normalized_abundance) %>% 
		spread(taxa, normalized_abundance) %>% 
		arrange(unique_id, day)

	surrogate_ts <- composite_ts %>% 
		group_by(unique_id) %>% 
		mutate(day = c(0,sample(day[day > 0]))) %>% 
		arrange(unique_id, day) %>% 
		ungroup
	
	best_E <- nonlinearity %>% 
		filter(taxa %in% c(i, otus)) %>% 
		select(taxa, embedding) %>% 
		unique

	print(paste0('Beginning interaction test of ', 
		paste(i, 'is driven by', paste(otus, collapse = ', ')), ' in ', treatment_subset))
# check mae for best theta
	test_theta <- lapply(1:500, function(blank){
		rnd_order <- data.frame(unique_id = sample(mouse_list, replace = T), 
			order = 1:length(mouse_list), stringsAsFactors = F)
		rnd_composite_ts <- inner_join(composite_ts, rnd_order, by = 'unique_id') %>% 
			arrange(order, day)

		theta_run <- block_lnlp(select(rnd_composite_ts, one_of(i, otus)),
			method = 's-map',
			num_neighbors = 0,
			theta = c(0, 1e-04, 3e-04, 0.001,
				0.003, 0.01, 0.03, 0.1,
                0.3, 0.5, 0.75, 1, 1.5,
                2, 3, 4, 6, 8),
			target_column = i,
			silent = T)
		return(theta_run)
	})
	test_theta <- do.call('rbind', test_theta)
#	test_theta %>% 
#		ggplot(aes(x = theta, y = mae)) +
#			geom_point(alpha = 0.01) +
#			stat_summary(fun.y = median, geom = 'line')
	best_theta <- test_theta %>% 
		group_by(theta) %>% 
		summarise(median_mae = median(mae, na.rm = T)) %>% 
		filter(median_mae == min(median_mae)) %>% 
		pull(theta)
# run smap with best theta
# partial derivaties of multivariate smap approximates interactions
	interaction_smap <- lapply(1:500, function(blank){
		rnd_order <- data.frame(unique_id = sample(mouse_list, replace = T), 
			order = 1:length(mouse_list), stringsAsFactors = F)
		rnd_composite_ts <- inner_join(composite_ts, rnd_order, by = 'unique_id') %>% 
			arrange(order, day)
		
		smap_res <-  block_lnlp(select(rnd_composite_ts, one_of(i, otus)),
			method = "s-map",
			num_neighbors = 0, 
			theta = best_theta,
			target_column = i,
			silent = T,
			save_smap_coefficients = T) # save S-map coefficients
		smap_coef_df <- smap_res$smap_coefficients[[1]]
		colnames(smap_coef_df) <- c(paste0('d', i, 'd', i), paste0('d', i, 'd', otus), 'intercept')
		return(bind_cols(smap_res$model_output[[1]],
			smap_coef_df,
			rnd_composite_ts))		
	})

	interaction_smap <- do.call('rbind', interaction_smap)

	interaction_smap %>% 
		ggplot(aes(x = obs, y = pred)) + 
			geom_point(alpha = 0.01) + 
			theme_bw()

	## Time series of fluctuating interaction strength

	# Plot all partial derivatives
	interaction_plot <- interaction_smap  %>% 
		gather(interaction, strength, one_of(paste0('d', i, 'd', otus))) %>% 
		group_by(unique_id, day) %>% 
		summarise(median_strength = median(strength)) %>% 
		ungroup %>% 
		mutate(time = 1:length(day)) %>% 
		ggplot(aes(x = time, y =median_strength)) + 
			geom_line() + 
			theme_bw()

	ggsave(filename = paste0(save_dir, '/interactions_w_', i, '.jpg'), 
		plot = interaction_plot, width = 7, height = 10, device = 'jpeg')
	
	write.table(interaction_smap, paste0(save_dir, '/interactions_w_', i, '.txt'), 
		quote = F, row.names = F)

}



# output file with taxa, interaction direction, interaction strength

