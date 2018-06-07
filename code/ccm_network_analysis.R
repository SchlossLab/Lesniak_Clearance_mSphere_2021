library(tidyverse)
library(statnet)
library(cowplot)

data_path <- "scratch/ccm_all/"   # path to the data
files <- dir(data_path, pattern = "ccm_raw_data*") # get file names
treatment_list <- unique(gsub('\\d{,2}.txt', '', files))

for(treatment in treatment_list){
	
	files <- dir(data_path, pattern = paste0(treatment, "*")) # get file names
	data <- files %>%
	  # read in all the files, appending the path before the filename
	  map(~ read.table(file.path(data_path, .), header = T, stringsAsFactors = F)) %>% 
	  reduce(rbind)

	data <- bind_rows(
		select(data, otu = otu1, strength = otu1_cause_otu2, p_value = pval_a_cause_b, affected_otu = otu2,
			prediction_slope = otu1_prediction_slope, p_slope = otu1_prediction_slope_p, E = E_A),
		select(data, otu = otu2, strength = otu2_cause_otu1, p_value = pval_b_cause_a, affected_otu = otu1,
			prediction_slope = otu2_prediction_slope, p_slope = otu2_prediction_slope_p, E = E_B))

	interaction_data <- data %>% 
		group_by(otu1, otu2) %>% 
		summarise(med_p_a_b = median(pval_a_cause_b),
			med_p_b_a = median(pval_b_cause_a),
			otu1_cause_otu2 = median(otu1_cause_otu2),
			otu2_cause_otu1 = median(otu2_cause_otu1)) %>% 
		gather(pval_comparison, pvalue, contains('med_p')) %>% 
		gather(interaction, strength, contains('cause')) %>% 
		#arrange(desc(otu1), desc(otu2)) %>%
		mutate(adj_strength = ifelse(pvalue > 0.05, 0, strength)) 

	interaction_data %>% 
		ggplot(aes(otu2, otu1)) + geom_tile(aes(fill = adj_strength)) + 
			scale_fill_gradient(low = 'white', high = 'blue')

test_data <- interaction_data %>% 
	filter(otu1 %in% c('Otu000001', 'Otu000003', 'Otu000004'), 
		otu2 %in% c('Otu000001', 'Otu000003', 'Otu000004'))

	num_nodes <- length(unique(data$otu1))
	network_matrix <- matrix()

}



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
	
significant_otus#