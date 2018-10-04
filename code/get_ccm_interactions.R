library(cowplot)
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
		median_rho_max_lobs > 0,
		p_min_v_max_lobs < 0.05,
		linear_corr_p < 0.05,
		ccm_null_p < 0.05) %>% 
	select(causal) %>% 
	separate(causal, sep = 'xmap', c('driven', 'driver'))
#i <- "C_difficile"
for(i in unique(xmap_otus$driven)){
	otus <- xmap_otus %>% 
		filter(driven == i) %>% 
		pull(driver)
	print(paste(i, 'is driven by', paste(otus, collapse = ', ')))

	composite_ts <- otu_df %>% 
		filter(taxa %in% c(i, otus), differenced == 'first') %>% 
		select(taxa, day, unique_id, normalized_abundance) %>% 
		spread(taxa, normalized_abundance) %>% 
		arrange(unique_id, day)		
	
	data_by_plot <- split(composite_ts, composite_ts$unique_id)
	segments_end <- cumsum(sapply(data_by_plot, NROW))
	segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
	segments <- cbind(segments_begin, segments_end) 
	pred_num <- round(0.2 * length(mouse_list))
	lib_segments <- segments[-tail(1:NROW(segments), pred_num), ]
	pred_segments <- segments[tail(1:NROW(segments), pred_num), ]

	best_E <- nonlinearity %>% 
		filter(taxa %in% c(i, otus)) %>% 
		select(taxa, embedding) %>% 
		unique
	univariate_ts <- bind_cols(composite_ts[, c('day', 'unique_id')], 
		make_block(composite_ts[,i],
			lib = segments))#,
			#max_lag = pull(filter(best_E, taxa == i), embedding))

	print(paste0('Beginning interaction test of ', 
		paste(i, 'is driven by', paste(otus, collapse = ', ')), ' in ', treatment_subset))
# check mae for best theta
	test_theta <- list()
	for(iter in 1:500){
		rndpred <- sample(1:length(mouse_list), pred_num)
		pred_order <- data.frame(
			unique_id = c(sample(mouse_list[-rndpred], replace = T), mouse_list[rndpred]),
			order = 1:NROW(mouse_list), 
			stringsAsFactors = F)
		rnd_composite_ts <- right_join(composite_ts, pred_order, by = 'unique_id')
		rnd_univariate_ts <- right_join(univariate_ts, pred_order, by = 'unique_id')

		uni_theta <- block_lnlp(select(univariate_ts, contains('col1')),
			method = 's-map',
			num_neighbors = 0,
			lib = lib_segments, pred = pred_segments,
			theta = c(0, 1e-04, 3e-04, 0.001,
				0.003, 0.01, 0.03, 0.1,
                0.3, 0.5, 0.75, 1, 1.5,
                2, 3, 4, 6, 8),
			target_column = 'col1',
			silent = T)
		theta_run <- block_lnlp(select(composite_ts, one_of(i, otus)),
			method = 's-map',
			num_neighbors = 0,
			lib = lib_segments, pred = pred_segments,
			theta = c(0, 1e-04, 3e-04, 0.001,
				0.003, 0.01, 0.03, 0.1,
                0.3, 0.5, 0.75, 1, 1.5,
                2, 3, 4, 6, 8),
			target_column = i,
			silent = T)
		test_theta[[iter]] <- bind_rows(
			data.frame(theta_run, model_type = 'multivariate', stringsAsFactors = F),
			data.frame(uni_theta, model_type = 'univariate', stringsAsFactors = F))
	}
	test_theta <- do.call('rbind', test_theta)
#	test_theta %>% 
#		ggplot(aes(x = theta, y = mae)) +
#			geom_point(alpha = 0.01) +
#			stat_summary(fun.y = median, geom = 'line')
	best_theta <- test_theta %>% 
		group_by(model_type, theta) %>% 
		summarise(median_mae = median(mae, na.rm = T)) %>% 
		filter(median_mae == min(median_mae)) 
	rm(test_theta)
# run smap with best theta
# partial derivaties of multivariate smap approximates interactions
	interaction_smap <- list()
	for(iter in 1:500){
		pred_order <- data.frame(unique_id = sample(mouse_list, replace = T),
			order = 1:NROW(mouse_list), stringsAsFactors = F)
		rnd_composite_ts <- right_join(composite_ts, pred_order, by = 'unique_id')
		rnd_univariate_ts <- right_join(univariate_ts, pred_order, by = 'unique_id')
		smap_multi <-  block_lnlp(select(rnd_composite_ts, one_of(i, otus)),
			lib = segments, pred = segments,
			method = "s-map",
			num_neighbors = 0, 
			theta = pull(filter(best_theta, model_type == 'multivariate'), theta),
			target_column = i,
			silent = T,
			save_smap_coefficients = T) # save S-map coefficients
		smap_uni <-  block_lnlp(select(univariate_ts, contains('col1')),
			lib = segments, pred = segments,
			method = "s-map",
			num_neighbors = 0, 
			theta = pull(filter(best_theta, model_type == 'univariate'), theta),
			target_column = 'col1',
			silent = T,
			save_smap_coefficients = T) # save S-map coefficients
		smap_coef_df <- smap_multi$smap_coefficients[[1]]
		colnames(smap_coef_df) <- c(paste0('d', i, '_d', i), paste0('d', i, '_d', otus), 'intercept')
		interaction_smap[[iter]] <- bind_rows(
			data.frame(bind_cols(smap_multi$model_output[[1]], smap_coef_df,
					right_join(composite_ts, pred_order, by = 'unique_id')), 
				embed = 'multi', stringsAsFactors = F),
			data.frame(smap_uni$model_output[[1]], embed = 'uni', stringsAsFactors = F))		
	}

	interaction_smap <- do.call('rbind', interaction_smap)

	write.table(interaction_smap, paste0(save_dir, '/interactions_w_', i, '.txt'), 
	quote = F, row.names = F)

	# comparison to univariate is limited to late time points 
	# so  limited interpretation
	#interaction_smap %>% 
	#	ggplot(aes(x = obs, y = pred, color = embed)) + 
	#		geom_point(alpha = .2) + 
	#		theme_bw() + 
	#	    coord_cartesian(ylim=c(-2, 2), xlim=c(-2, 2))
	#interaction_smap %>% 
	#	group_by(embed) %>% 
	#	summarise(mean = mean(obs, na.rm = T), median = median(obs, na.rm = T),
	#		pval = cor.test(.$obs, .$pred, method = 'spearman')$p.val,
	#		cor = cor.test(.$obs, .$pred, method = 'spearman')$estimate)

	## Time series of fluctuating interaction strength
	order <- expand.grid(unique_id = mouse_list, day = 0:10)  %>% 
		arrange(unique_id, day) %>% 
		mutate(epochs = 1:NROW(.)) %>% 
		right_join(expand.grid(epochs = 1:NROW(.), run = 1:500)) %>% 
		arrange(unique_id, run, day)
	# Plot all partial derivatives
	interaction_plot <- interaction_smap %>% 
		filter(embed == 'multi') %>% 
		right_join(order, by = c('day', 'unique_id')) %>% 
		gather(interaction, strength, one_of(paste0('d', i, '_d', otus))) %>% 
		mutate(interaction = gsub('tu0*', 'TU', interaction)) %>% 
		ggplot(aes(x = epochs, y = strength)) + 
			stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
				alpha = 0.2, fun.args =(conf.int = 0.5), aes(fill = interaction)) + 
			stat_summary(aes(color = interaction), fun.y = median, geom = 'line') + 
			theme_bw() +
			theme(legend.position="top", legend.title=element_blank()) +
			labs(title = paste0('Strength of effect of OTUs on ', gsub('tu0*', 'TU', i)))
	
	dynamics_plot <- shared_file %>% 
			select(cage, mouse, day, one_of(i, otus)) %>% 
			mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
			right_join(data.frame(unique_id = rep(order, each = 11), day = rep(0:10, length(order)))) %>% 
			mutate(time = 1:nrow(.)) %>% 
			gather(bacteria, counts, one_of(i, otus)) %>% 
				ggplot(aes(x = time, y = counts, color = as.factor(cage), group = interaction(cage, mouse))) + 
					geom_line(alpha = 0.4) + 
					geom_point() + 
					facet_grid(bacteria~., scales = 'free_y') +
					theme_bw() + 
					labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
						subtitle = 'Temporal Dynamics - Colored by mouse', 
						title = paste(i, 'is driven by', paste(otus, collapse = ', '))) + 
					theme_bw(base_size = 8) + 
					theme(legend.position = 'none', panel.grid.minor = element_blank())

	ggsave(filename = paste0(save_dir, '/interactions_w_', i, '.jpg'), 
		plot = plot_grid(dynamics_plot, interaction_plot, nrow = 2, align = 'v'),
		width = 14, height = 14, device = 'jpeg')
	

}



# output file with taxa, interaction direction, interaction strength

