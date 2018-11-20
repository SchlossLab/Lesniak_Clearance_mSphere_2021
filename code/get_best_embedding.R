library(rEDM)
library(tidyverse)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
input_file <- as.character(input_values[2])
save_dir <- paste0('data/process/ccm/')
print(paste0('Running set ', run_set))

#input_file <- 'data/process/ccm_otu_data.txt' 
abx_df <- read.table(input_file, header = T, stringsAsFactors = F)

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
abx_df <- abx_df %>% 
	filter(treatment == treatment_subset, differenced == 'first') 

print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset))

ifelse(!dir.exists(save_dir), dir.create(save_dir), print(paste0(save_dir, ' directory ready')))
ifelse(!dir.exists(paste0(save_dir, treatment_subset)), 
	dir.create(paste0(save_dir, treatment_subset)), 
	print(paste0(save_dir, treatment_subset, ' directory ready')))
ifelse(!dir.exists(paste0(save_dir, treatment_subset, '/embedding')), 
	dir.create(paste0(save_dir, treatment_subset, '/embedding')), 
	print(paste0(save_dir, treatment_subset, '/embedding directory ready')))

# remove otus that are present in less than 10 samples
taxa_list <- unique(abx_df$otu_feature)
sample_list <- unique(abx_df$unique_id) 

set.seed(seed)
# Choose random segments for prediction

taxa_nonlinearity_df <- list()
# test each otu for embedding and nonlinearity
# taxa <- 1
for(taxa in 1:length(taxa_list)){
	taxa_var <- taxa_list[[taxa]]

	composite_ts <- filter(abx_df, otu_feature == taxa_var)
	pred_num <- round(0.2 * length(sample_list))
	data_by_plot <- split(composite_ts, composite_ts$unique_id)
	segments_end <- cumsum(sapply(data_by_plot, NROW))
	segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
	segments <- cbind(segments_begin, segments_end) 
	lib_segments <- segments[-tail(1:NROW(segments), pred_num), ]
	pred_segments <- segments[tail(1:NROW(segments), pred_num), ]

	simplex_cat <- list()
	for(i in 1:500){ # when increasing seem to get repeats with treatment with 5 mice
		rndpred <- sample(1:length(sample_list), pred_num)
		pred_order <- data.frame(
			unique_id = c(sample(sample_list[-rndpred], replace = T), sample_list[rndpred]),
			order = 1:NROW(sample_list), 
			stringsAsFactors = F)
		rnd_composite_ts <- right_join(composite_ts, pred_order, by = 'unique_id')
		simplex_out <- simplex(data.frame(select(rnd_composite_ts, day, normalized_abundance)), 
			E = 2:6, lib = lib_segments, pred = pred_segments)#,

# inspect predictions of simplex		
#			stats_only = F)$model_output
#			simplex_out$mae
#			data.frame(select(composite_ts, day, normalized_abundance))[composite_pred[1]:composite_pred[2],]
		simplex_cat[[i]] <- cbind(simplex_out, run = i)
	}
	simplex_cat <- do.call('rbind', simplex_cat)

	simplex_median_mae <- simplex_cat %>% 
		group_by(E) %>% 
		summarise(median_mae = median(mae),
			lower_iqr_mae = quantile(mae)[2], 
			upper_iqr_mae = quantile(mae)[4])

	# set E by time series (not specific lib/pred split)
	best_E <- simplex_median_mae %>% 
		filter(median_mae == min(median_mae, na.rm = T)) %>% 
		pull(E)

	embedded_plot <- simplex_cat %>% 
		ggplot(aes(x = E, y = mae)) +
			geom_line(aes(group = run), alpha = 0.05) + 
			stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
				alpha = 0.2, fun.args =(conf.int = 0.5), 
				color = '#CC0000', linetype = 3, fill = NA) + 
			stat_summary(fun.y = median, geom = 'line', color = '#CC0000') + 
			geom_vline(xintercept = best_E, color = '#009999') + 
			labs(title = 'Simplex plot', subtitle = 'Selected embedding highlighted with vertical blue line\nMedian (red solid line) with IQR (red dashed lines)') + 
			theme_bw(base_size = 8)

	ggsave(filename = paste0(save_dir, treatment_subset, '/embedding/', taxa_var,  
		'_simplex_embedding_plot.jpg'),  
		plot = embedded_plot, 
		width = 7, height = 10, device = 'jpeg') 

	taxa_nonlinearity_df[[taxa]] <- data.frame(taxa = taxa_var, simplex_median_mae)


	print(paste('Completed ', taxa_var))
}

taxa_nonlinearity_df <- do.call('rbind', taxa_nonlinearity_df)

write.table(taxa_nonlinearity_df, paste0(save_dir, treatment_subset, '/simplex_embedding_first_differenced.txt'), 
	quote = F, row.names = F)

