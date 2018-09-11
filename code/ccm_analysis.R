library(rEDM)
library(tidyverse)
library(gtools)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
save_dir <- paste0('scratch/ccm_otu/')
print(paste0('Running set ', run_set))

ccm_otu_df <- 'data/process/ccm_otu_data.txt'
ccm_otu_df   <- read.table(ccm_otu_df, header = T, stringsAsFactors = F) %>% 
	mutate(otu_feature = gsub('_first', '', otu_feature))
#source('code/sum_otu_by_taxa.R')
#taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'
#shared_by_genus <- sum_otu_by_taxa(taxonomy_file = taxonomy_file, 
#	otu_df = shared_file, 
#	taxa_level = 'genus')

seed <- 062818
treatment_subset <- unique(ccm_otu_df$treatment)[run_set]
# use to test function
# most complete sample sets
#	treatment_subset <- 'cef_0.5_FALSE'
#	treatment_subset <- 'clinda_10_FALSE'

print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset))

ifelse(!dir.exists(paste0(save_dir, treatment_subset, '/ccm')), 
	dir.create(paste0(save_dir, treatment_subset, '/ccm')), 
	print(paste0(save_dir, treatment_subset, '/ccm directory ready')))

# check for embedding file
if(!file.exists(paste0(save_dir, treatment_subset, 
	'/smap_nonlinearity_first_differenced.txt'))){ 
	stop('Run run_smap.R first or find smap_nonlinearity_first_differenced.txt')
}
# load file with embeddings and nonlinearity tests for each otu
embedding_nonlinearity <- read.table(paste0(save_dir, treatment_subset, 
	'/smap_nonlinearity_first_differenced.txt'), header = T, stringsAsFactors = F) %>% 
	mutate(taxa = gsub('_first', '', taxa))

abx_df <- ccm_otu_df %>% 
	filter(treatment == treatment_subset) %>% 
	filter(otu_feature %in% unique(embedding_nonlinearity$taxa))

rm(ccm_otu_df)

# create list of mice and taxa
taxa_list <- best_embedding$taxa
mouse_list <- unique(abx_df$unique_id)

# create a list of all combinations of taxa
otu_combinations <- apply(combinations(length(taxa_list), 2, repeats=TRUE), 1, list)

set.seed(seed)

n_samples <- nrow(filter(abx_df, otu_feature == 'C_difficile_first'))
lib_sizes <- c(seq(5, n_samples, by = 5))

run_ccm <- function(otu, input_df, treatment_subset, taxa_list){
	current_otu1 <- taxa_list[ otu[[1]][1] ]
	current_otu2 <- taxa_list[ otu[[1]][2] ]
	
	composite_ts <- abx_df %>% 
		select(unique_id, day, normalized_abundance, otu_feature) %>% 
		filter(otu_feature %in% c(current_otu1, current_otu2)) %>% 
		spread(otu_feature, normalized_abundance)
	surrogate_ts <- composite_ts %>% 
		group_by(unique_id) %>% 
		mutate(day = c(0,sample(day[day > 0]))) %>% 
		arrange(unique_id, day) %>% 
		ungroup
	data_by_plot <- split(composite_ts, composite_ts$unique_id)
	segments_end <- cumsum(sapply(data_by_plot, NROW))
	segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
	segments <- cbind(segments_begin, segments_end) 
	
	best_E <- embedding_nonlinearity %>% 
		filter(taxa %in% c(current_otu1, current_otu2)) %>% 
		select(taxa, embedding) %>% 
		unique

	print(paste0('Beginning ', current_otu1, ' and ', current_otu2, ' in from ', treatment_subset))

	ccm_run_results <- lapply(1:500, function(i){
		rnd_order <- sample(1:NROW(lib_segments), replace = T)
		rnd_segments <- segments[rnd_order, ]
		mice_order <- paste(mouse_list[rnd_order], collapse = '--')
		
		#Run the CCM test
		# Does B "cause" A?
		otu1_xmap_otu2 <- ccm(composite_ts, lib = rnd_segments, pred = rnd_segments, 
			lib_column = current_otu1, target_column = current_otu2, 
			E = pull(filter(best_E, taxa %in% current_otu1), embedding), 
			lib_sizes = lib_sizes, silent = TRUE)
		# Does A "cause" B?
		otu2_xmap_otu1 <- ccm(composite_ts, lib = rnd_segments, pred = rnd_segments, 
			lib_column = current_otu2, target_column = current_otu1, 
			E = pull(filter(best_E, taxa %in% current_otu2), embedding), 
    		lib_sizes = lib_sizes, silent = TRUE)
		# null model with surrogate data
		otu1_xmap_otu2_surr <- ccm(surrogate_ts, lib = rnd_segments, pred = rnd_segments, 
			lib_column = current_otu1, target_column = current_otu2, 
			E = pull(filter(best_E, taxa %in% current_otu1), embedding), 
			lib_sizes = lib_sizes, silent = TRUE)
		otu2_xmap_otu1_surr <- ccm(surrogate_ts, lib = rnd_segments, pred = rnd_segments, 
			lib_column = current_otu2, target_column = current_otu1, 
			E = pull(filter(best_E, taxa %in% current_otu2), embedding), 
    		lib_sizes = lib_sizes, silent = TRUE)

		return(rbind(
			data.frame(rbind(otu1_xmap_otu2, otu2_xmap_otu1), data = 'real'),
			data.frame(rbind(otu1_xmap_otu2_surr, otu2_xmap_otu1_surr), data = 'surrogate'))
			)
		})

	ccm_data <- do.call('rbind', ccm_run_results) %>% 
		mutate_if(is.factor, as.character)

	# plot the ability of otu to predict the other otu
	# note causal direction is opposite of xmap
	# so if otu1 xmaps to otu2, then otu2 causes otu1
	# this is because otu2 has left an impression on otu1 

	default_rho <- cor.test(pull(composite_ts, current_otu1), pull(composite_ts, current_otu2),
		method = "pearson", use = 'complete.obs')

	ccm_plot <- ccm_data$rho %>% 
		mutate(causal = paste0(lib_column, ' xmap ', target_column)) %>% 
		ggplot(aes(x = lib_size, y = rho)) +
			stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
				alpha = 0.2, fun.args =(conf.int = 0.5), aes(fill = data)) + 
			stat_summary(aes(color = data), fun.y = median, geom = 'line') + 
			scale_color_manual(values = c('#CC0000', '#555555'), limits = c('real', 'surrogate')) +
			scale_fill_manual(values = c('#CC0000', '#555555'), limits = c('real', 'surrogate')) + 
			guides(color = FALSE, fill=guide_legend(title=NULL)) + 
			labs(x = 'theta', y = 'delta mae', title = 'S-map Analysis - Test feature nonlinearity', 
				subtitle = 'Surrogate data is randomly permuted time indices\nMedian (solid line) with IQR (dotted lines) and 5/95th percentile (shaded area)\n(First differenced, treatment = Antibiotic_Dose_RecoveryBeforeChallenge)') + 
			facet_grid(causal~., scale = 'free_y') + 
			geom_hline(yintercept = default_rho$estimate, linetype = 3) + 
			theme_bw(base_size = 8) + theme(legend.position = c(0.1,0.9))

	ggsave(filename = paste0(save_dir, treatment_subset, '/ccm/', current_otu1, 
		'_', current_otu2, 'ccm.jpg'), 
		plot = ccm_plot, width = 7, height = 10, device = 'jpeg')
	
	# test for convergence
		# test for significant monotonic increasing trend with rho(L) using kendall's tau test
	kendalls_test_otu1 <- cor(
		pull(filter(ccm_data, lib_column %in% current_otu1), rho),
		pull(filter(ccm_data, lib_column %in% current_otu1), lib_size),
		method = 'kendall')
	kendalls_test_otu2 <- cor(
		pull(filter(ccm_data, lib_column %in% current_otu2), rho),
		pull(filter(ccm_data, lib_column %in% current_otu2), lib_size),
		method = 'kendall')

		# test for significant improvement  in rho(L) by Fisher's delta rho Z test 
		# (test for difference min/max library)
	fishers_z_test <- 

	# test for significance of cross map
	linear_corr_test <- 

	null_test <- 

	taxa_nonlinearity_df <- rbind(taxa_nonlinearity_df, 
		data.frame(taxa = taxa_var, embedding = best_E, nonlinear_output, real_v_surrogate_nonlinearity))

	print(paste('Completed ', taxa_var))
}

write.table(taxa_nonlinearity_df, paste0(save_dir, treatment_subset, '/smap_nonlinearity_first_differenced.txt'), 
	quote = F, row.names = F)



	lobs_test_range <- round(length(unique(ccm_plot_df$lobs))*0.1)

	ccm_data <- ccm_plot_df %>% 
		group_by(driver_otu) %>% 
		mutate(time_point = case_when(lobs %in% head(sort(unique(lobs)), lobs_test_range) ~ 'initial',
			lobs %in% tail(sort(unique(lobs)), lobs_test_range) ~ 'end',
			T ~ 'middle')) %>% 
		filter(time_point != 'middle') %>% 	
		summarise(ccm_p_value = wilcox.test(rho~time_point, alternative = 'greater')$p.value) %>%
		full_join(ccm_data, by = c('driver_otu'))

	title <- ggdraw() + 
	  draw_label(paste0(treatment_subset, ' with ', current_otu1, ' and ', current_otu2,
	  	'\n(Data is first differenced)\n(treatment = Antibiotic_Dose_Allow recovery before C difficile Challenge)'),
		fontface = 'bold')
	if(min(ccm_data$ccm_p_value) < 0.05/10^5){ # ~100,000 comparisons across all treatments
		ggsave(filename = paste0(save_dir, treatment_subset, '/sig_ccm_', current_otu1, 
				'_', current_otu2, '_first_differenced.jpg'),
			plot = plot_grid(title, plot_grid(plot_grid(lagged_dynamics_plot, dynamics_plot, embedding_dim_plot, prediction_step_plot), 
				CCM_plot, align = 'v', ncol = 1, labels = 'AUTO'),  ncol = 1, rel_heights = c(0.1, 1)),
			width = 7, height = 10, device = 'jpeg')
		} else {
		ggsave(filename = paste0(save_dir, treatment_subset, '/ccm_', current_otu1, 
				'_', current_otu2, '_first_differenced.jpg'),
			plot = plot_grid(title, plot_grid(plot_grid(lagged_dynamics_plot, dynamics_plot, embedding_dim_plot, prediction_step_plot), 
				CCM_plot, align = 'v', ncol = 1, labels = 'AUTO'),  ncol = 1, rel_heights = c(0.1, 1)),
			width = 7, height = 10, device = 'jpeg')
	}
	
	print(paste0('Completed ', current_otu1, ' and ', current_otu2,  ' from ', treatment_subset))
	return(ccm_data)
}

print(paste0('Beginning Treatment Set - ', treatment_subset, ' (Antibiotic, Dosage, Delay Challenge with C difficile)'))


print('Beginning CCM on 1st differenced data')

output <- map_df(otu_combinations, ~ run_ccm(., input_df = data.frame(abx_df), 
	treatment_subset = treatment_subset, taxa_list = taxa_list))
write.table(output, paste0(save_dir, treatment_subset, '/ccm_by_genus_', treatment_subset, '_first_differenced_', seed, 'seed.txt'), 
	quote = F, row.names = F)

print(paste0('Completed treatment set - ', treatment_subset))