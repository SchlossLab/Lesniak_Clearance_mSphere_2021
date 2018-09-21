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

}

# for all interactions that are significant
# import embeddings

# check mae for best theta
# run smap with best theta
# partial derivaties of multivariate smap approximates interactions

# output file with taxa, interaction direction, interaction strength

