library(rEDM)
library(tidyverse)
library(gtools)
library(methods)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
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
save_dir <- paste0('scratch/ccm_otu/', treatment_subset, '/ccm')
ifelse(!dir.exists(save_dir), 
	dir.create(save_dir), 
	print(paste0(save_dir, ' directory ready')))
ifelse(!dir.exists(paste0(save_dir, '/temp')), 
	dir.create(paste0(save_dir, '/temp')), 
	print(paste0(save_dir, '/temp directory ready')))

# check for embedding file
if(!file.exists(paste0(save_dir, '/../smap_nonlinearity_first_differenced.txt'))){ 
	stop('Run run_smap.R first or find smap_nonlinearity_first_differenced.txt')
}
# load file with embeddings and nonlinearity tests for each otu
embedding_nonlinearity <- read.table(paste0(save_dir, 
	'/../smap_nonlinearity_first_differenced.txt'), header = T, stringsAsFactors = F) %>% 
	mutate(taxa = gsub('_first', '', taxa)) %>% 
	filter(p_linear_v_nonlinear < 0.05, p_real_v_surrogate < 0.05, data == 'real')

abx_df <- ccm_otu_df %>% 
	filter(treatment == treatment_subset) %>% 
	filter(otu_feature %in% unique(embedding_nonlinearity$taxa))

rm(ccm_otu_df)

# create list of mice and taxa
taxa_list <- unique(embedding_nonlinearity$taxa)
mouse_list <- unique(abx_df$unique_id)

# create a list of all combinations of taxa
otu_combinations <- apply(combinations(length(taxa_list), 2, repeats=TRUE), 1, list)


n_samples <- nrow(filter(abx_df, otu_feature == unique(embedding_nonlinearity$taxa)[1]))
lib_sizes <- c(seq(5, n_samples, by = 5))

run_ccm <- function(otu, input_df, treatment_subset, taxa_list){
	set.seed(seed)
	# current_otu1 <- taxa_list[ 1 ]
	# current_otu2 <- taxa_list[ 1 ]
	current_otu1 <- taxa_list[ otu[[1]][1] ]
	current_otu2 <- taxa_list[ otu[[1]][2] ]
	if(!file.exists(paste0(save_dir, '/temp/ccm_by_otu_', current_otu1, '_', 
		current_otu2, '_', treatment_subset, '_first_differenced.txt'))){ 
		# input_df <- abx_df
		composite_ts <- input_df %>% 
			select(unique_id, day, normalized_abundance, otu_feature) %>% 
			filter(otu_feature %in% c(current_otu1, current_otu2)) %>% 
			spread(otu_feature, normalized_abundance)

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
			surrogate_ts <- composite_ts %>% 
				group_by(unique_id) %>% 
				mutate(day = c(0,sample(day[day > 0]))) %>% 
				arrange(unique_id, day) %>% 
				ungroup
			pred_order <- data.frame(unique_id = sample(mouse_list, replace = T),
				order = 1:NROW(mouse_list), stringsAsFactors = F)
			rnd_composite_ts <- right_join(composite_ts, pred_order, by = 'unique_id')
			rnd_surrogate_ts <- right_join(surrogate_ts, pred_order, by = 'unique_id')
			#Run the CCM test
			# Does B "cause" A?
			otu1_xmap_otu2 <- ccm(rnd_composite_ts, lib = segments, pred = segments, 
				lib_column = current_otu1, target_column = current_otu2, 
				E = pull(filter(best_E, taxa %in% current_otu1), embedding), 
				lib_sizes = lib_sizes, silent = TRUE)
			# Does A "cause" B?
			otu2_xmap_otu1 <- ccm(rnd_composite_ts, lib = segments, pred = segments, 
				lib_column = current_otu2, target_column = current_otu1, 
				E = pull(filter(best_E, taxa %in% current_otu2), embedding), 
	    		lib_sizes = lib_sizes, silent = TRUE)
			# null model with surrogate data
			otu1_xmap_otu2_surr <- ccm(surrogate_ts, lib = segments, pred = segments, 
				lib_column = current_otu1, target_column = current_otu2, 
				E = pull(filter(best_E, taxa %in% current_otu1), embedding), 
				lib_sizes = lib_sizes, silent = TRUE)
			otu2_xmap_otu1_surr <- ccm(surrogate_ts, lib = segments, pred = segments, 
				lib_column = current_otu2, target_column = current_otu1, 
				E = pull(filter(best_E, taxa %in% current_otu2), embedding), 
	    		lib_sizes = lib_sizes, silent = TRUE)

			return(rbind(
				data.frame(rbind(otu1_xmap_otu2, otu2_xmap_otu1), data = 'real'),
				data.frame(rbind(otu1_xmap_otu2_surr, otu2_xmap_otu1_surr), data = 'surrogate'))
				)
			})

		ccm_data <- do.call('rbind', ccm_run_results) %>% 
			mutate_if(is.factor, as.character) %>% 
			mutate(causal = paste0(lib_column, 'xmap', target_column))

		# plot the ability of otu to predict the other otu
		# note causal direction is opposite of xmap
		# so if otu1 xmaps to otu2, then otu2 causes otu1
		# this is because otu2 has left an impression on otu1 

		default_rho <- cor.test(pull(composite_ts, current_otu1), pull(composite_ts, current_otu2),
			method = "spearman", use = 'complete.obs')

		ccm_plot <- ccm_data %>% 
			ggplot(aes(x = lib_size, y = rho)) +
				stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
					alpha = 0.2, fun.args =(conf.int = 0.5), aes(fill = data)) + 
				stat_summary(aes(color = data), fun.y = median, geom = 'line') + 
				scale_color_manual(values = c('#CC0000', '#555555'), limits = c('real', 'surrogate')) +
				scale_fill_manual(values = c('#CC0000', '#555555'), limits = c('real', 'surrogate')) + 
				guides(color = FALSE, fill=guide_legend(title=NULL)) + 
				labs(x = 'libray size', y = 'rho', title = 'Convergent Cross Map Analysis - Test feature interaction', 
					subtitle = 'Surrogate data is randomly permuted time indices\nMedian (solid line) with IQR (dotted lines) and 5/95th percentile (shaded area)\n(First differenced, treatment = Antibiotic_Dose_RecoveryBeforeChallenge)') + 
				facet_grid(causal~., scale = 'free_y') + 
				geom_hline(yintercept = default_rho$estimate, linetype = 3) + 
				theme_bw(base_size = 8) + theme(legend.position = c(0.1,0.9))

		ggsave(filename = paste0(save_dir, '/', current_otu1, '_', current_otu2, '_ccm.jpg'), 
			plot = ccm_plot, width = 7, height = 10, device = 'jpeg')
		
		# test for convergence
			# test for significant monotonic increasing trend with rho(L) using spearman test
			# test for significant improvement in rho(L) by wilcox test (difference min/max library)
		print('Running ccm convergence test')
		ccm_convergence_test <- full_join( 
			group_by(ccm_data, causal, data) %>% 
				summarise(ccm_trend_rho = cor.test(rho, lib_size, 
					method = 'spearman', use = 'complete.obs')$estimate,
					ccm_trend_p = cor.test(rho, lib_size, 
					method = 'spearman', use = 'complete.obs')$p.value),
			filter(ccm_data, lib_size %in% c(min(lib_size), max(lib_size))) %>%
				group_by(causal, data) %>% 
				summarise(median_rho_min_lobs = median(rho[lib_size == min(lib_size)], na.rm = T),
					median_rho_max_lobs = median(rho[lib_size == max(lib_size)], na.rm = T),
					p_min_v_max_lobs = wilcox.test(rho ~ lib_size,
						alternative = 'less')$p.value),
			by = c('causal','data')) %>% 
			ungroup

		# test for significance of cross map
		print('Running ccm significance test')
		ccm_significance_test <- full_join(
			group_by(ccm_data, causal, data) %>% 
				filter(lib_size > quantile(lib_size)['75%']) %>% 
				summarise(linear_corr_p = my.t.test.p.value(rho, 
						mu = default_rho$estimate, alternative = 'greater')) %>% 
				ungroup,
			group_by(ccm_data, causal) %>% 
				filter(lib_size > quantile(lib_size)['75%']) %>% 
				summarise(ccm_null_p = wilcox.test(rho ~ data, 
						alternative = 'greater')$p.value) %>% 
				ungroup,
			by = 'causal')
		print('Saving ccm results')
		write.table(full_join(ccm_convergence_test, ccm_significance_test, by = c("causal", "data")), 
			paste0(save_dir, '/temp/ccm_by_otu_', 
				current_otu1, '_', current_otu2, '_', treatment_subset, '_first_differenced.txt'), 
		quote = F, row.names = F)

		print(paste0('Completed ', current_otu1, ' and ', current_otu2, ' in from ', 
			treatment_subset))
	} else {
		print(paste0(current_otu1, ' and ', current_otu2, ' in from ', 
			treatment_subset, ' already complete'))
	}

}

# create function to generate NA if t test cannot calculate value
my.t.test.p.value <- function(...) {
	obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

print(paste0('Beginning Treatment Set - ', treatment_subset, ' (Antibiotic, Dosage, Delay Challenge with C difficile)'))

print('Beginning CCM on 1st differenced data')

for(x in otu_combinations){ 
	run_ccm(x, input_df = data.frame(abx_df), 
		treatment_subset = treatment_subset, taxa_list = taxa_list)
}

output_temp <- list.files(paste0(save_dir, '/temp/'))
ccm_cat_df <- do.call("rbind", 
	lapply(paste0(save_dir, '/temp/', output_temp), read.table, header = TRUE)) 
write.table(ccm_cat_df,
		paste0(save_dir, '/../ccm_by_otu_', treatment_subset, '_first_differenced.txt'), 
	quote = F, row.names = F)

unlink(paste0(save_dir, '/temp/'), recursive = T)

print(paste0('Completed treatment set - ', treatment_subset))