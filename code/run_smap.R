library(rEDM)
library(tidyverse)
library(cowplot)
library(gtools)
#library(viridis)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
save_dir <- paste0('scratch/ccm_otu/')
print(paste0('Running set ', run_set))

ccm_otu_df <- 'data/process/ccm_otu_data.txt'
abx_df <- read.table(ccm_otu_df, header = T, stringsAsFactors = F)
#source('code/sum_otu_by_taxa.R')
#taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'
#shared_by_genus <- sum_otu_by_taxa(taxonomy_file = taxonomy_file, 
#	otu_df = shared_file, 
#	taxa_level = 'genus')

seed <- 062818
treatment_subset <- unique(abx_df$treatment)[run_set]
# use to test function
# most complete sample sets
#	treatment_subset <- 'cef_0.5_FALSE'
#	treatment_subset <- 'clinda_10_FALSE'
ifelse(!dir.exists(paste0(save_dir, treatment_subset, '/nonlinearity')), 
	dir.create(paste0(save_dir, treatment_subset, '/nonlinearity')), 
	print(paste0(save_dir, treatment_subset, '/nonlinearity directory ready')))

abx_df <- abx_df %>% 
	filter(treatment == treatment_subset) 
# check for embedding file
if(!file.exists(paste0(save_dir, treatment_subset, 
	'/simplex_embedding_first_differenced.txt'))){ 
	stop('Run get_best_embedding.R first or find simplex_embedding_first_differenced.txt')
}
# load file with embeddings for each otu
embedding <- read.table(paste0(save_dir, treatment_subset, 
	'/simplex_embedding_first_differenced.txt'), header = T, stringsAsFactors = F)

# currently selecting lowest median mae
# could also try to select lowest E which is not significantly
# different from lowest median
best_embedding <- embedding %>% 
	group_by(taxa) %>% 
	filter(median_mae == min(median_mae)) %>% 
	filter(E == min(E)) %>% 
	ungroup	

print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset))

# remove otus that are present in less than 10 samples
taxa_list <- best_embedding$taxa
mouse_list <- unique(abx_df$unique_id) 
theta <- seq(0, 2, 0.1)

# create a list of all combinations of taxa
otu_combinations <- apply(combinations(length(taxa_list), 2, repeats=TRUE), 1, list)

set.seed(seed)

taxa_nonlinearity_df <- c()
# test each otu for nonlinearity
# taxa_var <- taxa_list[[30]]
for(taxa_var in taxa_list){
	s_map_cat <- c()
	for(i in 1:1000){
		best_E <- filter(best_embedding, taxa == taxa_var) %>% 
			pull(E)
		composite_ts <- filter(abx_df, otu_feature == taxa_var)
		surrogate_ts <- composite_ts %>% 
			group_by(unique_id) %>% 
			mutate(day = c(0,sample(day[day > 0]))) %>% 
			arrange(unique_id, day) %>% 
			ungroup
		data_by_plot <- split(composite_ts, composite_ts$unique_id)
		segments_end <- cumsum(sapply(data_by_plot, NROW))
		segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
		segments <- cbind(segments_begin, segments_end) 
		rndpred <- sample(1:NROW(segments), 1)
		lib_segments <- segments[-rndpred, ]
		rndlib <- sample(1:NROW(lib_segments), replace = T)
		composite_lib <- lib_segments[rndlib, ]
		composite_pred <- segments[rndpred, ]
		smap_out <- s_map(data.frame(select(composite_ts, day, normalized_abundance)), 
			E = best_E, lib = composite_lib, pred = composite_pred, theta = theta)
		surrogate_smap <- s_map(data.frame(select(surrogate_ts, day, normalized_abundance)), 
			E = best_E, lib = composite_lib, pred = composite_pred, theta = theta)
		s_map_cat <- rbind(s_map_cat, 
			rbind(cbind(smap_out, data = 'real', run = i),
				cbind(surrogate_smap, data = 'surrogate', run = i)))
	}

	s_map_cat <- s_map_cat %>% 
		select(run, mae, theta, data) %>% 
		arrange(data, run, theta) %>% 
		group_by(data, run) %>% 
		mutate(linear_mae = head(mae, 1)) %>% 
		group_by(data, run, theta) %>% 
		mutate(delta_mae = mae - linear_mae) %>% 
		ungroup
		#group_by(data, theta) %>% 
		#summarise(delta_mae = max(delta_mae, na.rm = T)) %>% 
	smap_plot <- s_map_cat %>% 
		ggplot(aes(x = theta, y = delta_mae)) +
			#geom_line(alpha = 0.2)
			stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
				alpha = 0.2, fun.args =(conf.int = 0.5), linetype = 3,
				aes(color = data, fill = NA)) + 
			stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
				alpha = 0.2, fun.args =(conf.int = 0.9), aes(fill = data)) + 
			stat_summary(aes(color = data), fun.y = median, geom = 'line') + 
			scale_color_manual(values = c('#CC0000', '#555555'), limits = c('real', 'surrogate')) +
			scale_fill_manual(values = c('#CC0000', '#555555'), limits = c('real', 'surrogate')) + 
			guides(color = FALSE, fill=guide_legend(title=NULL)) + 
			labs(x = 'theta', y = 'delta mae', title = 'S-map Analysis - Test feature nonlinearity', 
				subtitle = 'Surrogate data is randomly permuted time indices\nMedian (solid line) with IQR (dotted lines) and 5/95th percentile (shaded area)') + 
			theme_bw(base_size = 8) + theme(legend.position = c(0.1,0.9))

	nonlinear_output <-  full_join(
			group_by(s_map_cat, data) %>% 
			filter(theta == max(theta)) %>% 
			summarise(median_delta_mae = median(delta_mae)), 
	# test each for decrease in mae (increase in prediciton with increased theta)
			group_by(s_map_cat, data) %>% 
			filter(theta %in% c(min(theta), max(theta))) %>%
			summarise(p_linear_v_nonlinear = wilcox.test(mae ~ theta,
				alternative = 'greater')$p.value),
		by = 'data')

	# test if mae is lower in real vs surrogate data
	real_v_surrogate_nonlinearity <- s_map_cat %>% 
		filter(theta == max(theta)) %>% 
		summarise(p_real_v_surrogate = wilcox.test(mae ~ data,
				alternative = 'less')$p.value)

	title <- ggdraw() + 
	  draw_label(paste0(treatment_subset, ' with ', taxa_var,
	  	'\n(First differenced, treatment = Antibiotic_Dose_RecoveryBeforeChallenge)'),
		fontface = 'bold')
	ggsave(filename = paste0(save_dir, treatment_subset, '/nonlinearity/', taxa_var, 
		'_simplex_smap.jpg'), 
		plot = plot_grid(title, smap_plot,  ncol = 1, rel_heights = c(0.1, 1)),
			width = 7, height = 10, device = 'jpeg')
	taxa_nonlinearity_df <- rbind(taxa_nonlinearity_df, 
		data.frame(taxa = taxa_var, embedding = best_E, delta_mae))

	print(paste('Completed ', taxa_var))
}

nonlinear_taxa <- taxa_nonlinearity_df %>% 
	filter(median_delta_mae > 0, data == 'real') %>% 
	mutate(taxa_diff = taxa) %>% 
	separate(taxa, c('otu', 'differenced')) %>% 
	mutate(diff = case_when(differenced == 'raw' ~ 0,
		differenced == 'first' ~ 1,
		differenced == 'second' ~ 2)) %>% 
	group_by(otu) %>% 
	filter(diff == min(diff))