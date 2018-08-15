library(multispatialCCM)
library(tidyverse)
library(cowplot)
library(gtools)
library(viridis)
library(forecast)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
save_dir <- paste0('scratch/ccm_otu/')
print(paste0('Running set ', run_set))

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(group, cage, mouse, day, CFU, cdiff, abx, dose, delayed) %>% 
	unite(treatment, abx, dose, delayed) %>% 
	filter(cdiff == T, day >= 0, treatment != 'none_NA_FALSE')
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F)
#source('code/sum_otu_by_taxa.R')
#taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'
#shared_by_genus <- sum_otu_by_taxa(taxonomy_file = taxonomy_file, 
#	otu_df = shared_file, 
#	taxa_level = 'genus')

seed <- 062818
treatment_subset <- unique(meta_file$treatment)[run_set]
# use to test function
# most complete sample sets
#	treatment_subset <- 'cef_0.5_FALSE'
#	treatment_subset <- 'clinda_10_FALSE'

print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset))

ifelse(!dir.exists(save_dir), dir.create(save_dir), print(paste0(save_dir, ' directory ready')))
ifelse(!dir.exists(paste0(save_dir, treatment_subset)), 
	dir.create(paste0(save_dir, treatment_subset)), 
	print(paste0(save_dir, treatment_subset, ' directory ready')))

abx_df <- meta_file %>% 
	filter(treatment == treatment_subset) %>%
	mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
	inner_join(shared_file, by = c('group' = "Group")) %>% 
	select(-group, -numOtus) %>% 
	rename(C_difficile = CFU)
		
# remove otus that are present in less than 10 samples
abx_df <- select(abx_df, day, C_difficile, which(apply(abx_df > 1, 2, sum) > 10 )) 
taxa_list <- colnames(select(abx_df, -day, -cage, -mouse, -treatment, -unique_id))

normalize <- function(x, ...) {
    (x - mean(x, ...))/sd(x, ...)
}

abx_df <- abx_df %>% 
	# replace all 0s with random value between 0 and 1
	gather(taxa, abundance, one_of(taxa_list)) %>% 
	mutate(abundance = ifelse(abundance == 0, 
		sample(100, sum(abundance == 0), replace = T)/100, abundance)) %>% 
	spread(taxa, abundance) %>% 
	# add missing days
	full_join(data.frame(
		unique_id = rep(unique(abx_df$unique_id), each = 11),
		day = rep(seq(0,10), length(unique(abx_df$unique_id))),
		stringsAsFactors = F), by = c('unique_id', 'day')) %>% 
	# impute values for missing days and normalize by mouse
	arrange(unique_id, day) %>% 
	gather(bacteria, raw_abundance, C_difficile, contains('Otu00')) %>% 
	group_by(unique_id, bacteria) %>% 
	nest() %>% 
	mutate(abundance = map(data, ~normalize(na.interp(.$raw_abundance)))) %>% 
	unnest(data, abundance) %>% 
	select(unique_id, day, bacteria, abundance) 

# create a list of all combinations of taxa
otu_combinations <- apply(combinations(length(taxa_list), 2, repeats=TRUE), 1, list)

# use to test function
#input_df <- abx_df
#otu <- otu_combinations[[2]]

setup_df_for_mccm <- function(input_df){
	# reorder mice
	sample_mice <- sample(unique(input_df$unique_id), replace = F) 
	n_mice <- length(sample_mice)
	output_df <- data.frame(unique_id = sample_mice, sample = 1:n_mice, stringsAsFactors = F) %>% 
		inner_join(input_df, by = 'unique_id') %>% 
		arrange(bacteria, sample, day) %>%  
		spread(bacteria, abundance)
	# set day 0 to NA to separate data by mouse for ccm (for 1st differenced, day 0 == NA)
	output_df[which(output_df$day == 0), ] <- NA
	return(list(abundance_df = output_df, mice_order = sample_mice))
}

run_ccm <- function(otu, input_df, treatment_subset, taxa_list){
	current_otu1 <- taxa_list[ otu[[1]][1] ]
	current_otu2 <- taxa_list[ otu[[1]][2] ]

 	input_df <- input_df %>%  
	# create a 1st differenced dataframe
		arrange(unique_id, bacteria, day) %>% 
		group_by(unique_id, bacteria) %>% 
		mutate_at(vars(abundance) , funs(. - lag(.))) %>% 
		ungroup

	set.seed(seed)
	print(paste0('Beginning ', current_otu1, ' and ', current_otu2, ' in from ', treatment_subset))

	ccm_run_results <- lapply(1:500, function(i){
		ccm_df <- setup_df_for_mccm(input_df)
		Accm <- pull(ccm_df$abundance_df, current_otu1)
		Bccm <- pull(ccm_df$abundance_df, current_otu2)
		mice_order <- paste(ccm_df$mice_order, collapse = '--')
Accm[is.na(Accm)] <- 0
test_nonlinearity(Accm, e=3)

data(e120_invnit16)
str(e120_invnit16)
e120_invnit16 %>% 
	ggplot(aes(x = Year, y = noh020tot, group = as.factor(Plot), color = as.factor(Plot))) + 
		geom_line()

simplex_output <- simplex(Accm, E = 2, tp = 1:10)
plot(simplex_output$tp, simplex_output$rho, type = "l", xlab = "Time to Prediction (tp)", 
	ylab = "Forecast Skill (rho)")
smap_output <- s_map(Accm, E = 2)
plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", 
	ylab = "Forecast Skill (rho)")

		#Maximum E to test - one less than number of observations per sample
		# ideal to be at minimum E or lower dim, prevent overfitting by selecting lower dim with moderate pred power
		maxE<- 7 #length(unique(abx_df$day)) - 2 # one less for separating NAs and one less sample
		#Matrix for storing output
		Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c(current_otu1, current_otu2)
		#Loop over potential E values and calculate predictive ability
		#of each process for its own dynamics
		for(E in 2:maxE) {
		#Uses defaults of looking forward one prediction step (predstep)
		#And using time lag intervals of one time step (tau)
		Emat[E-1,1]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
		Emat[E-1,2]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
		}
		#maximum E 
		# ideal to be at minimum E or lower dim, prevent overfitting by selecting lower dim with moderate pred power
		maxEmat <- Emat/c(2:maxE)
		E_A<- c(2:maxE)[which(maxEmat[,1] == max(maxEmat[,1], na.rm =T))]
		E_B<- c(2:maxE)[which(maxEmat[,2] == max(maxEmat[,2], na.rm =T))]
		
		# create df for embedding dim plot
		embedding_dim_df <- data.frame(cbind(Emat, E = c(2:maxE))) %>% 
				gather(bacteria, rho, -E) %>% 
				left_join(data.frame(bacteria = c(current_otu1, current_otu2), 
					Selected_E = c(E_A, E_B),
					run = i, stringsAsFactors = F ), by = 'bacteria')

		#Check data for nonlinear signal that is not dominated by noise
		#Checks whether predictive ability of processes declines with
		#increasing time distance
		#See manuscript and R code for details
		signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
		predsteplist=1:10)
		signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
		predsteplist=1:10)

		pred_plot_df <- rbind(data.frame(signal_A_out$predatout, bacteria = current_otu1, run = i),
				data.frame(signal_B_out$predatout, bacteria = current_otu2, run = i))
		
		#Run the CCM test
		#E_A and E_B are the embedding dimensions for A and B.
		#tau is the length of time steps used (default is 1)
		#iterations is the number of bootsrap iterations (default 100)
		# 100 iterations is sufficient to reduce the range in performance
		# Does A "cause" B?
		CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=10)
		# Does B "cause" A?
		CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=10)
		ccm_plot_df <- rbind(data.frame(driver_otu = current_otu1,
					driven_otu = current_otu2,
					run = i,
					gather(data.frame(CCM_boot_A$FULLinfo, lobs = CCM_boot_A$Lobs),
						iter, rho, -lobs)),
				data.frame(driver_otu = current_otu2,
					driven_otu = current_otu1,
					run = i,
					gather(data.frame(CCM_boot_B$FULLinfo, lobs = CCM_boot_B$Lobs),
						iter, rho, -lobs)))
		#Test for significant causal signal
		#See R function for details
		CCM_significance_test<-ccmtest(CCM_boot_A, CCM_boot_B)
		ccm_data <- rbind(data.frame(
				run = i,  
				driver_otu = current_otu1,
				driven_otu = current_otu2,
				driver_predicts_driven = max(CCM_boot_A$rho), 
				# pval = 1 - ( # of iterations final rho > last rho / total iterations)
				portion_noncausal = CCM_significance_test[['pval_a_cause_b']], 
				E = E_A,
				# slope and p determined using summary(lm())
				self_prediction_slope = signal_A_out$rho_pre_slope['Estimate'],
				slope_p = signal_A_out$rho_pre_slope['Pr(>|t|)'],
				treatment = treatment_subset, mice = mice_order, stringsAsFactors = F), 
			data.frame(
				run = i,  
				driver_otu = current_otu2,
				driven_otu = current_otu1,
				driver_predicts_driven = max(CCM_boot_B$rho),
				# pval = 1 - ( # of iterations final rho > last rho / total iterations)
				portion_noncausal = CCM_significance_test[['pval_b_cause_a']],
				E = E_B,
				# slope and p determined using summary(lm())
				self_prediction_slope = signal_B_out$rho_pre_slope['Estimate'],
				slope_p = signal_B_out$rho_pre_slope['Pr(>|t|)'],
				treatment = treatment_subset, mice = mice_order, stringsAsFactors = F)) %>% 
			mutate(non_linear = ifelse(slope_p < ( 0.05/length(run) ), T, F)) %>% 
			separate(treatment, c('abx', 'dose', 'delayed_infection'), sep = '_')
		return(list(embed = embedding_dim_df, pred = pred_plot_df, ccm = ccm_plot_df, data = ccm_data))
	})
	print('CCM completed, generating plots now')

	embedding_dim_df <- do.call('rbind', lapply(ccm_run_results, '[[', 'embed')) %>% 
		mutate_if(is.factor, as.character)
	pred_plot_df <- do.call('rbind', lapply(ccm_run_results, '[[', 'pred')) %>% 
		mutate_if(is.factor, as.character)
	ccm_plot_df <- do.call('rbind', lapply(ccm_run_results, '[[', 'ccm')) %>% 
		mutate_if(is.factor, as.character)
	ccm_data <- do.call('rbind', lapply(ccm_run_results, '[[', 'data')) %>% 
		mutate_if(is.factor, as.character)

	# plot each time point agasint the previous day
	lagged_dynamics_plot <- input_df %>% 
			mutate(t1 = lag(abundance)) %>% 
			filter(!is.na(t1), !is.na(abundance)) %>% 
			ggplot(aes(x = abundance, y = t1, color = day)) + 
				geom_point() + facet_wrap(~bacteria, scales = 'free') +  
				theme_bw(base_size = 8) + 
				labs(x= expression('t'[0]), y = expression('t'[1]))
	# plot temporal dynamics of otus
	dynamics_plot <- meta_file %>% 
		filter(treatment == treatment_subset) %>%
		inner_join(shared_file, by = c('group' = "Group")) %>%
		rename(C_difficile = CFU) %>% 
		select(cage, mouse, day, one_of(current_otu1, current_otu2)) %>% 
		gather(bacteria, counts, one_of(current_otu1, current_otu2)) %>% 
			ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), group = interaction(cage, mouse))) + 
				geom_line(alpha = 0.4) + 
				geom_point() + 
				facet_grid(bacteria~., scales = 'free_y') +
				theme_bw() + 
				labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
					title = 'Temporal Dynamics', subtitle = 'Colored by mouse') + 
				scale_x_continuous(breaks=seq(0,10, 1)) + 
				theme_bw(base_size = 8) + 
				theme(legend.position = 'none', panel.grid.minor = element_blank())
	# plot embedding dimension of each otu/sample with the indicated used value for E
	embedding_dim_plot <- embedding_dim_df %>% 
		ggplot(aes(x = E, y = rho, color = bacteria, group = interaction(bacteria, run))) + 
			geom_line(alpha = 0.1) + 
			geom_vline(aes(xintercept = Selected_E, color = bacteria), 
				linetype = 'dashed', size = 0.5, alpha = 0.3, show.legend = FALSE) +
			labs(x = 'E', y = 'Pearson correlation coefficient (rho)', title = 'Embedding Dimension Selection',
				subtitle = 'Dimension of highest predictive power') + 
			theme_bw(base_size = 8) + 
			theme(legend.position = c(0.2, 0.2), legend.title=element_blank(), 
				legend.background=element_blank(), panel.grid.minor = element_blank()) + 
			scale_x_continuous(breaks = seq(2, max(embedding_dim_df$E), 1)) + 
			guides(colour = guide_legend(override.aes = list(alpha = 1)))
	# plot prediction over time, to determine if prediction decays with time (indicative of non-linearity)
	prediction_step_plot <- pred_plot_df %>% 
		left_join(select(ccm_data, non_linear, run, driver_otu), 
			by = c('bacteria' = 'driver_otu', 'run')) %>% 
		ggplot(aes(x = predstep, y = rho, color = non_linear, group = interaction(bacteria, run))) + 
			geom_line(alpha = 0.3) + facet_grid(bacteria~., scales = 'free_y') + 
			labs(x = 'Prediction Steps', y = 'Pearson correlation coefficient (rho)', 
				title = 'Predictive Power', color = 'Nonlinear') + 
			theme_bw(base_size = 8) + 
			theme(legend.position = c(0.9, 0.9), legend.background=element_blank(),
				panel.grid.minor = element_blank()) + 
			scale_x_continuous(breaks = seq(1, 10, 1)) +
			scale_color_manual(values = c('gray', 'red'))
	# plot the ability of otu to predict the other otu
	CCM_plot <- ccm_plot_df %>% 
		left_join(select(ccm_data, non_linear, run, driver_otu), 
				by = c('driver_otu', 'run')) %>% 
		mutate(causal = paste0(driver_otu, ' causes ', driven_otu)) %>% 
		ggplot(aes(x = lobs, y = rho)) +
				geom_hex(bins = 100) + 
				scale_fill_viridis() + 
				facet_grid(.~causal) + 
				theme(legend.position = 'none') + 
				labs(x = 'L', y = 'Pearson correlation coefficient (rho)', color = '', fill = '') + 
				theme_bw() + 
				theme(legend.position='none') + 
				guides(colour = guide_legend(override.aes = list(alpha = 1)))
	
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