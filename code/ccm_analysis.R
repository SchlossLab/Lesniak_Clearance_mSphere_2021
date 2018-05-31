library(multispatialCCM)
library(tidyverse)
library(cowplot)

#		test workflow with short simulated RPS data
#		setup to work with diff abx combinations
#		check into ientification of otus with 0 abundance
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

seed <- commandArgs(TRUE)
print(paste0('Using ', seed, ' as seed'))
set.seed(seed)

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(group, cage, mouse, day, CFU, cdiff, abx, dose, delayed) %>% 
	unite(treatment, abx, dose, delayed) %>% 
	filter(cdiff == T, day >= 0, treatment != 'none_NA_FALSE')
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_file <- read.table(shared_file, sep = '\t', header = T)
output <- c()

print(paste0('Beginning seed ', seed))

ifelse(!dir.exists('scratch/ccm'), dir.create('scratch/ccm'), print('ccm/ directory ready'))

for(treatment_subset in unique(meta_file$treatment)){
	print(paste0('Beginning Treatment Set - ', treatment_subset, ' (Antibiotic, Dosage, Delay Challenge with C difficile)'))
	# treatment_subset <- 'amp_0.5_TRUE' # test for missing day 0 (missing 1)
	# treatment_subset <- 'amp_0.5_FALSE' # test for missing day 0 (missing 2)
	# treatment_subset <- 'cef_0.1_FALSE' # test for missing day 0 (none missing)

	abx_df <- meta_file %>% 
		filter(treatment == treatment_subset) %>%
		mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
		inner_join(select(shared_file, -label, -numOtus),
			by = c('group' = "Group"))
	
	abx_df <- select(abx_df, day, CFU, which(apply(abx_df > 1, 2, sum) > 10 ))

	missing_day_0 <- summarise(group_by(abx_df, unique_id), first_day = min(day)) %>% 
		filter(first_day != 0) %>% 
		pull(unique_id)

	if(length(missing_day_0) > 0){
		abx_df <- abx_df %>% 
			bind_rows(data.frame(unique_id = missing_day_0, day = 0)) %>% 
			mutate(random_order = as.numeric(factor(unique_id, levels = 
				sample(unique(unique_id), length(unique(unique_id)), replace = F)))) %>% 
			arrange(random_order, day) %>%  
			select(day, CFU, contains('Otu'))
		} else {
		abx_df <- abx_df %>% 
			mutate(random_order = as.numeric(factor(unique_id, levels = 
				sample(unique(unique_id), length(unique(unique_id)), replace = F)))) %>% 
			arrange(random_order, day) %>%  
			select(day, CFU, contains('Otu'))
		}

	NA_list <- which(abx_df$day == 0)
	zero_subset <- abx_df == 0 & !is.na(abx_df)
	abx_df[zero_subset] <-  sample(100,sum(zero_subset), replace = T)/100
	abx_df[NA_list, ] <- NA

	# create folder for treatment set
	ifelse(!dir.exists(file.path('scratch/ccm', treatment_subset)), 
		dir.create(file.path('scratch/ccm', treatment_subset)), FALSE)

	for(i in which(grepl('Otu', colnames(abx_df)))){
		Accm<-abx_df$CFU
		Bccm<-abx_df[,i]
		current_otu <- colnames(abx_df)[i]
		print(paste0('Beginning ', current_otu, ' in from ', treatment_subset))
		lagged_dynamics_plot <- data.frame(day = abx_df$day,
			t1 = abx_df[,i],
			t0 = c(abx_df[-1,i], NA))  %>% 
			filter(!is.na(t0), !is.na(t1)) %>% 
			ggplot(aes(x = t0, y = t1, color = day)) + 
				geom_point() + 
				labs(title = paste(current_otu)) + 
				theme_bw(base_size = 8)
		#Maximum E to test - one less than number of observations per sample
		# ideal to be at minimum E or lower dim, prevent overfitting by selecting lower dim with moderate pred power
		maxE<- 8 #length(unique(abx_df$day)) - 2 # one less for separating NAs and one less sample
		#Matrix for storing output
		Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("C_difficile", current_otu)
		#Loop over potential E values and calculate predictive ability
		#of each process for its own dynamics
		for(E in 2:maxE) {
		#Uses defaults of looking forward one prediction step (predstep)
		#And using time lag intervals of one time step (tau)
		Emat[E-1,"C_difficile"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
		Emat[E-1,current_otu]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
		}
		#maximum E 
		# ideal to be at minimum E or lower dim, prevent overfitting by selecting lower dim with moderate pred power
		maxEmat <- Emat/c(2:maxE)
		E_A<-c(2:maxE)[which(maxEmat[,1] == max(maxEmat[,1], na.rm =T))]
		E_B<-c(2:maxE)[which(maxEmat[,2] == max(maxEmat[,2], na.rm =T))]
		embedding_dim_plot <- data.frame(cbind(Emat, E = c(2:maxE))) %>% 
			gather(bacteria, rho, -E) %>% 
			left_join(data.frame(bacteria = c('C_difficile', current_otu), Selected_E = c(E_A, E_B))) %>% 
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
		prediction_step_plot <- rbind(data.frame(signal_A_out$predatout, bacteria = 'C_difficile'),
			data.frame(signal_B_out$predatout, bacteria = current_otu)) %>% 
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
		#Note - increase iterations to 100 for consistant results
		CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=100)
		# Does B "cause" A?
		CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=100)
		#Test for significant causal signal
		#See R function for details
		CCM_significance_test<-ccmtest(CCM_boot_A, CCM_boot_B)
		current_ccm <- data.frame(t(CCM_significance_test), 
			cdiff_cause_otu = max(CCM_boot_A$rho), 
			otu_cause_cdiff = max(CCM_boot_B$rho),
			otu = colnames(abx_df)[i],
			E_A = E_A,
			E_B = E_B,
			treatment = treatment_subset) %>% 
			separate(treatment, c('abx', 'dose', 'delayed_infection'), sep = '_')

		if(current_ccm$pval_b_cause_a <= 0.1){
			causal_otu <- as.character(current_ccm$otu)
			CCM_plot <- rbind(data.frame(causal = paste0('Cdiff_causes_', causal_otu), 
				lobs = CCM_boot_A$Lobs,
				rho = CCM_boot_A$rho,
				stdev_min = CCM_boot_A$rho - CCM_boot_A$sdevrho,
				stdev_max = CCM_boot_A$rho + CCM_boot_A$sdevrho),
			data.frame(causal = paste0(causal_otu, '_causes_Cdiff'), 
				lobs = CCM_boot_B$Lobs,
				rho = CCM_boot_B$rho,
				stdev_min = CCM_boot_B$rho - CCM_boot_B$sdevrho,
				stdev_max = CCM_boot_B$rho + CCM_boot_B$sdevrho )) %>% 
			#gather(level, value, rho, stdev_min, stdev_max) %>% 
			ggplot(aes(x = lobs)) + 
				geom_line(aes(y = rho, color = causal)) + 
				geom_ribbon(aes(ymin = stdev_min, ymax = stdev_max, fill = causal), alpha = 0.3) + 
				labs(x = 'L', y = 'Pearson correlation coefficient (rho)', color = '', fill = '') + 
				theme_bw() + 
				theme(legend.position="top", legend.direction="horizontal")

			dynamics_plot <- meta_file %>% 
				filter(treatment == treatment_subset) %>%
				inner_join(select(shared_file, -label, -numOtus),
					by = c('group' = "Group")) %>%
				select(cage, mouse, day, CFU, one_of(causal_otu)) %>% 
				gather(bacteria, counts, CFU, one_of(causal_otu)) %>% 
					ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), group = interaction(cage, mouse))) + 
						geom_line() + 
						facet_grid(bacteria~., scales = 'free_y') +
						theme_bw() + 
						labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
							title = 'Temporal Dynamics', subtitle = 'Colored by mouse') + 
						scale_x_continuous(breaks=seq(0,10, 1)) + 
						theme_bw(base_size = 8) + 
						theme(legend.position = 'none')

			ggsave(paste0('scratch/ccm/', treatment_subset, '/ccm_cdiff_caused_by_', causal_otu, '_seed', seed, '.jpg'),
				plot_grid(plot_grid(lagged_dynamics_plot, dynamics_plot, embedding_dim_plot, prediction_step_plot), 
					CCM_plot, align = 'v', ncol = 1, labels = 'AUTO'))
		}
		print(paste0('Completed ', current_otu, ' from ', treatment_subset))
		output <- rbind(output, current_ccm)
	}
	print(paste0('Completed treatment set - ', treatment_subset))
}
print(paste0('Completed seed ', seed))
write.table(output, paste0('scratch/ccm/ccm_raw_data_seed', seed, '.txt'), quote = F, row.names = F)


ccm_output <- read.table('scratch/ccm/ccm_output.txt', header = T)
significant_output <- ccm_output %>% 
	filter(pval_b_cause_a < 0.05 | pval_a_cause_b < 0.05) %>% 
	unite(treatment, abx, dose, delayed_infection, otu) %>% 
	filter(treatment == 'amp_0_5_Otu000011')
	group_by(treatment) %>% 
	summarise(cdiff_causes_otu = mean(cdiff_cause_otu),
		p_cdiff_causes_otu = mean(pval_a_cause_b),
		otu_causes_cdiff = mean(otu_cause_cdiff),
		p_otu_causes_cdiff = mean(pval_b_cause_a),
		number_of_repeats_detected = n()) 

significant_otus <- rbind(
	mutate(
		select(significant_output, 
			treatment,
			interaction = cdiff_causes_otu, 
			p_value = p_cdiff_causes_otu,
			number_of_repeats_detected),
		direction = c('cdiff_causes_otu')),
	mutate(
		select(significant_output, 
			treatment,
			interaction = otu_causes_cdiff, 
			p_value = p_otu_causes_cdiff,
			number_of_repeats_detected),
		direction = c('otu_causes_cdiff'))
	) %>% 
	filter(p_value < 0.05, number_of_repeats_detected > 5) %>% 
	separate(treatment, c('abx', 'dose', 'delayed_infection')
	
significant_otus