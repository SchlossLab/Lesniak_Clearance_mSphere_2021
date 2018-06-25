library(multispatialCCM)
library(tidyverse)
library(cowplot)

#		test workflow with short simulated RPS data
#		setup to work with diff abx combinations
#		check into identification of otus with 0 abundance
#		randomize mouse order in combining samples
# what to do with mice with incomplete samples (missing days)?
#	cipro_10_FALSE - day 5 missing 1;
#	vanco_0.625_FALSE - day 3 and 4 missing 1;
#	vanco_0.3_FALSE - day 2 missing 1, day 3 missing 3;
#	vanco_0.1_FALSE - day 3 and 9 missing 1, day 2 missing 2;
#	strep_0.5_false - day 5 missing 4;
#	metro_1_false - day 9 and 10 missing 1, day 2 missing 3;
#	clinda_10_false - day 4 missing 1;
#	amp_0.5_false - day 4 missing 1, day 0 missing 2, day 9 missing 4;
#	cef_0.3_false - day 4 and 7 missing 1;
#	strep_5_false - day 2 missing 1;
#	amp_0.5_true - day 0 and 4 missing 1
input_values <- commandArgs(TRUE)
run_set <- input_values[1]
set_E <- input_values[2]
save_dir <- paste0('scratch/ccm/')
print(paste0('Running set ', run_set))

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(group, cage, mouse, day, CFU, cdiff, abx, dose, delayed) %>% 
	unite(treatment, abx, dose, delayed) %>% 
	filter(cdiff == T, day >= 0, treatment != 'none_NA_FALSE')
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F)
source('code/sum_otu_by_taxa.R')
taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'
shared_by_genus <- sum_otu_by_taxa(taxonomy_file = taxonomy_file, 
	otu_df = shared_file, 
	taxa_level = 'genus')

seed_treatment <- expand.grid(seed = 1:10, treatment = unique(meta_file$treatment))[run_set, ]
seed <- seed_treatment$seed
treatment_subset <- as.character(seed_treatment$treatment)

print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset, ' using seed ', seed, ', and setting E to ', set_E))

set.seed(seed)

print(paste0('Beginning seed ', seed))

ifelse(!dir.exists(save_dir), dir.create(save_dir), print(paste0(save_dir, ' directory ready')))
ifelse(!dir.exists(paste0(save_dir, treatment_subset)), 
	dir.create(paste0(save_dir, treatment_subset)), 
	print(paste0(paste0(save_dir, treatment_subset), ' directory ready')))

# use to test function
#otu <- list(c(1), c(15))
#input_df <- abx_df
#input_df <- abx_df_1diff

run_ccm <- function(otu, input_df, treatment_subset, data_diff){
	Accm<- select(input_df, -day)[ , otu[[1]] ]
	Bccm<- select(input_df, -day)[ , otu[[2]] ]
#	Accm_diff1 <- select(abx_df_1diff, -day)[ , otu[[1]] ]
#	Bccm_diff1 <- select(abx_df_1diff, -day)[ , otu[[2]] ]
	current_otu1 <- colnames(select(input_df, -day))[ otu[[1]] ]
	current_otu2 <- colnames(select(input_df, -day))[ otu[[2]] ]
	print(paste0('Beginning ', current_otu1, ' and ', current_otu2, ' in from ', treatment_subset))

	lagged_dynamics_plot <- input_df %>% 
			select(day, one_of(current_otu1, current_otu2)) %>% 
			gather(otu, t0, -day) %>% 
			mutate(t1 = lag(t0)) %>% 
			filter(!is.na(t1), !is.na(t0)) %>% 
			ggplot(aes(x = t0, y = t1, color = day)) + 
				geom_point() + facet_wrap(~otu, scales = 'free') +  
				theme_bw(base_size = 8)

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

	embedding_dim_plot <- data.frame(cbind(Emat, E = c(2:maxE))) %>% 
		gather(bacteria, rho, -E) %>% 
		left_join(data.frame(bacteria = c(current_otu1, current_otu2), Selected_E = c(E_A, E_B))) %>% 
		ggplot(aes(x = E, y = rho, color = bacteria)) + 
			geom_line() + 
			geom_vline(aes(xintercept = Selected_E, color = bacteria), 
				linetype = 'dashed', size = 0.5, show.legend = FALSE) +
			labs(x = 'E', y = 'Pearson correlation coefficient (rho)', title = 'Embedding Dimension Selection',
				subtitle = 'Dimension of highest predictive power') + 
			theme_bw(base_size = 8) + 
			theme(legend.position = c(0.8, 0.8), legend.title=element_blank(), 
				legend.background=element_blank()) + 
			scale_x_continuous(breaks = seq(2, maxE, 1))

	#Check data for nonlinear signal that is not dominated by noise
	#Checks whether predictive ability of processes declines with
	#increasing time distance
	#See manuscript and R code for details
	signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
	predsteplist=1:10)
	signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
	predsteplist=1:10)

	prediction_step_plot <- rbind(data.frame(signal_A_out$predatout, bacteria = current_otu1),
		data.frame(signal_B_out$predatout, bacteria = current_otu2)) %>% 
		ggplot(aes(x = predstep, y = rho, color = bacteria)) + 
			geom_line() + 
			labs(x = 'Prediction Steps', y = 'Pearson correlation coefficient (rho)', 
				title = 'Predictive Power') + 
			theme_bw(base_size = 8) + 
			theme(legend.position = c(0.8, 0.8), legend.title=element_blank(), 
				legend.background=element_blank()) + 
			scale_x_continuous(breaks = seq(1, 10, 1))

	#Run the CCM test
	#E_A and E_B are the embedding dimensions for A and B.
	#tau is the length of time steps used (default is 1)
	#iterations is the number of bootsrap iterations (default 100)
	# Does A "cause" B?
	#Note - increase iterations to 1000 for consistant results
	CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=1000)
	# Does B "cause" A?
	CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=1000)
	#Test for significant causal signal
	#See R function for details
	CCM_significance_test<-ccmtest(CCM_boot_A, CCM_boot_B)
	current_ccm <- data.frame( 
		otu1 = current_otu1,
		otu2 = current_otu2,
		otu1_cause_otu2 = max(CCM_boot_A$rho), 
		otu2_cause_otu1 = max(CCM_boot_B$rho),
		pval_otu1_cause_otu2 = CCM_significance_test[['pval_a_cause_b']],
		pval_otu2_cause_otu1 = CCM_significance_test[['pval_b_cause_a']],
		E_A = E_A,
		E_B = E_B,
		otu1_prediction_slope = paste(signal_A_out$rho_pre_slope['Estimate']),
		otu1_prediction_slope_p = paste(signal_A_out$rho_pre_slope['Pr(>|t|)']),
		otu2_prediction_slope = paste(signal_B_out$rho_pre_slope['Estimate']),
		otu2_prediction_slope_p = paste(signal_B_out$rho_pre_slope['Pr(>|t|)']),
		treatment = treatment_subset) %>% 
		separate(treatment, c('abx', 'dose', 'delayed_infection'), sep = '_')

#	if(current_ccm$pval_b_cause_a <= 0.1){		
	causal_otu <- c(current_otu1, current_otu2)
	CCM_plot <- rbind(data.frame(causal = paste0(current_otu1, '_causes_', current_otu2), 
		lobs = CCM_boot_A$Lobs,
		rho = CCM_boot_A$rho,
		stdev_min = CCM_boot_A$rho - CCM_boot_A$sdevrho,
		stdev_max = CCM_boot_A$rho + CCM_boot_A$sdevrho),
	data.frame(causal = paste0(current_otu2, '_causes_', current_otu1), 
		lobs = CCM_boot_B$Lobs,
		rho = CCM_boot_B$rho,
		stdev_min = CCM_boot_B$rho - CCM_boot_B$sdevrho,
		stdev_max = CCM_boot_B$rho + CCM_boot_B$sdevrho)) %>% 
		#gather(level, value, rho, stdev_min, stdev_max) %>% 
		ggplot(aes(x = lobs)) + 
			#geom_line(aes(y = rho, color = causal)) + 
			geom_ribbon(aes(ymin = stdev_min, ymax = stdev_max, fill = causal), alpha = 0.2) + 
			geom_point(aes(y = rho, color = causal), alpha = 0.4) + 
			labs(x = 'L', y = 'Pearson correlation coefficient (rho)', color = '', fill = '') + 
			theme_bw() + 
			theme(legend.position="top", legend.direction="horizontal")

	dynamics_plot <- meta_file %>% 
		filter(treatment == treatment_subset) %>%
		inner_join(shared_by_genus, by = c('group' = "Group")) %>%
		select(cage, mouse, day, C_difficile = CFU, one_of(current_otu1, current_otu2)) %>% 
		gather(bacteria, counts, one_of(current_otu1, current_otu2)) %>% 
			ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), group = interaction(cage, mouse))) + 
				geom_line() + 
				facet_grid(bacteria~., scales = 'free_y') +
				theme_bw() + 
				labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
					title = 'Temporal Dynamics', subtitle = 'Colored by mouse') + 
				scale_x_continuous(breaks=seq(0,10, 1)) + 
				theme_bw(base_size = 8) + 
				theme(legend.position = 'none')

	title <- ggdraw() + 
	  draw_label(paste0(treatment_subset, ' with ', current_otu1, ' and ', causal_otu,
	  	'\n(Data is ', data_diff, ', using seed', seed,
	  	')\n(treatment = Antibiotic_Dose_Allow recovery before C difficile Challenge)'),
		fontface = 'bold')

	ggsave(filename = paste0(save_dir, treatment_subset, '/ccm_', current_otu1, 
			'_', causal_otu, '_', data_diff, '_seed', seed, '.jpg'),
		plot = plot_grid(title, plot_grid(plot_grid(lagged_dynamics_plot, dynamics_plot, embedding_dim_plot, prediction_step_plot), 
			CCM_plot, align = 'v', ncol = 1, labels = 'AUTO'),  ncol = 1, rel_heights = c(0.1, 1)),
		width = 7, height = 10, device = 'jpeg')
	
	print(paste0('Completed ', current_otu1, ' and ', current_otu2,  ' from ', treatment_subset))
	return(current_ccm)
}

#run_each_treatment <- function(treatment_subset){
	print(paste0('Beginning Treatment Set - ', treatment_subset, ' (Antibiotic, Dosage, Delay Challenge with C difficile)'))
	# treatment_subset <- 'amp_0.5_TRUE' # test for missing day 0 (missing 1)
	# treatment_subset <- 'amp_0.5_FALSE' # test for missing day 0 (missing 2)
	# treatment_subset <- 'cef_0.1_FALSE' # test for missing day 0 (none missing)

	abx_df <- meta_file %>% 
		filter(treatment == treatment_subset) %>%
		mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
		inner_join(shared_by_genus, by = c('group' = "Group"))
	
# remove otus that are present in less than 10 samples
	abx_df <- select(abx_df, day, CFU, which(apply(abx_df > 1, 2, sum) > 10 )) 
	taxa_list <- colnames(select(abx_df, -day, -CFU, -cage, -mouse, -treatment, -unique_id))
# create a 1st differenced dataframe
	abx_df_1diff <- abx_df %>% 
		arrange(unique_id, day) %>% 
		group_by(unique_id) %>% 
		mutate_at(vars(c('CFU', taxa_list)) , funs(. - lag(.))) %>% 
		ungroup

# find which mice are missing data for day 0
	missing_day_0 <- summarise(group_by(abx_df, unique_id), first_day = min(day)) %>% 
		filter(first_day != 0) %>% 
		pull(unique_id)

# add day 0 back to those missing and randomize the order of the mice
	randomize_order <- function(input_df){
		if(length(missing_day_0) > 0){
			output_df <- input_df abx_df_1diff %>% 
				bind_rows(data.frame(unique_id = missing_day_0, day = 0, stringsAsFactors = F)) %>% 
				mutate(random_order = as.numeric(factor(unique_id, levels = 
					sample(unique(unique_id), length(unique(unique_id)), replace = F)))) %>% 
				arrange(random_order, day) %>%  
				select(day, C_difficile = CFU, one_of(taxa_list))#contains('Otu'))

		} else {
			output_df <- input_df %>% 
				mutate(random_order = as.numeric(factor(unique_id, levels = 
					sample(unique(unique_id), length(unique(unique_id)), replace = F)))) %>% 
				arrange(random_order, day) %>%  
				select(day, C_difficile = CFU, one_of(taxa_list))#contains('Otu'))
		}
		NA_list <- which(output_df$day == 0)
		zero_subset <- output_df == 0 & !is.na(output_df)
		output_df[zero_subset] <-  sample(100,sum(zero_subset), replace = T)/100
		output_df[NA_list, ] <- NA
		#output_df <- select(output_df, -day)
		return(output_df)
	}

abx_df <- randomize_order(abx_df)
abx_df_1diff <- randomize_order(abx_df_1diff)

otu_combinations <- cross2(1:ncol(abx_df), 1:ncol(abx_df))

output <- map_df(otu_combinations, ~ run_ccm(., input_df = abx_df, treatment_subset = treatment_subset, data_diff = "not_differenced"))
write.table(output, paste0(save_dir, treatment_subset, '/ccm_by_genus_', treatment_subset, '_not_differenced_', seed, 'seed.txt'), 
	quote = F, row.names = F)

output <- map_df(otu_combinations, ~ run_ccm(., input_df = abx_df_1diff, treatment_subset = treatment_subset, data_diff = 'first_differenced'))
write.table(output, paste0(save_dir, treatment_subset, '/ccm_by_genus_', treatment_subset, '_first_differenced_', seed, 'seed.txt'), 
	quote = F, row.names = F)

print(paste0('Completed treatment set - ', treatment_subset))

print(paste0('Completed seed ', seed))


#
#
#ccm_output <- read.table('scratch/ccm/ccm_output.txt', header = T)
#significant_output <- ccm_output %>% 
#	filter(pval_b_cause_a < 0.05 | pval_a_cause_b < 0.05) %>% 
#	unite(treatment, abx, dose, delayed_infection, otu) %>% 
#	filter(treatment == 'amp_0_5_Otu000011')
#	group_by(treatment) %>% 
#	summarise(cdiff_causes_otu = mean(cdiff_cause_otu),
#		p_cdiff_causes_otu = mean(pval_a_cause_b),
#		otu_causes_cdiff = mean(otu_cause_cdiff),
#		p_otu_causes_cdiff = mean(pval_b_cause_a),
#		number_of_repeats_detected = n()) 
#
#significant_otus <- rbind(
#	mutate(
#		select(significant_output, 
#			treatment,
#			interaction = cdiff_causes_otu, 
#			p_value = p_cdiff_causes_otu,
#			number_of_repeats_detected),
#		direction = c('cdiff_causes_otu')),
#	mutate(
#		select(significant_output, 
#			treatment,
#			interaction = otu_causes_cdiff, 
#			p_value = p_otu_causes_cdiff,
#			number_of_repeats_detected),
#		direction = c('otu_causes_cdiff'))
#	) %>% 
#	filter(p_value < 0.05, number_of_repeats_detected > 5) %>% 
#	separate(treatment, c('abx', 'dose', 'delayed_infection')
#	
#significant_otus#