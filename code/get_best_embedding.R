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
mouse_list <- unique(abx_df$unique_id) 

# create a list of all combinations of taxa
otu_combinations <- apply(combinations(length(taxa_list), 2, repeats=TRUE), 1, list)

set.seed(seed)
# Choose random segments for prediction

taxa_nonlinearity_df <- c()
# test each otu for embedding and nonlinearity
for(taxa_var in taxa_list){
	simplex_cat <- c()
	for(i in 1:500){ # when increasing seem to get repeats with treatment with 5 mice
		composite_ts <- filter(abx_df, otu_feature == taxa_var)
		data_by_plot <- split(composite_ts, composite_ts$unique_id)
		segments_end <- cumsum(sapply(data_by_plot, NROW))
		segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
		segments <- cbind(segments_begin, segments_end) 
		rndpred <- sample(1:NROW(segments), 1)
		lib_segments <- segments[-rndpred, ]
		rndlib <- sample(1:NROW(lib_segments), replace = T)
		composite_lib <- lib_segments[rndlib, ]
		composite_pred <- segments[rndpred, ]
		simplex_out <- simplex(data.frame(select(composite_ts, day, abundance)), 
			E = 2:6, lib = composite_lib, pred = composite_pred)
		simplex_cat <- rbind(simplex_cat, cbind(simplex_out, run = i))
	}
	names(simplex_out) <- taxa_var

	# set E by time series (not specific lib/pred split)
	best_E <- simplex_cat %>% 
		group_by(E) %>% 
		summarise(median_rho = median(rho, na.rm = T)) %>% 
		#mutate(adj_rho = median_rho / E) %>% 
		filter(median_rho == max(median_rho, na.rm = T)) %>% 
		pull(E)

	embedded_plot <- simplex_cat %>% 
		ggplot(aes(x = E, y = rho)) +
		#geom_smooth()
		geom_line(aes(group = run), alpha = 0.1) + 
		geom_point(alpha = 0.1) + 
		geom_vline(xintercept = best_E, color = 'red') +
		labs(title = 'Simplex plot', subtitle = 'Selected embedding highlighted with red line') + 
		theme_bw(base_size = 8)

	ggsave(filename = paste0(save_dir, treatment_subset, '/embedding/', taxa_var,  
		'_simplex_embedding_plot.jpg'),  
		plot = embedded_plot, 
		width = 7, height = 10, device = 'jpeg') 

	taxa_nonlinearity_df <- rbind(taxa_nonlinearity_df, 
		data.frame(taxa = taxa_var, embedding = best_E))


	print(paste('Completed ', taxa_var))
}

write.table(taxa_nonlinearity_df, paste0(save_dir, treatment_subset, '/simplex_embedding_first_differenced.txt'), 
	quote = F, row.names = F)

