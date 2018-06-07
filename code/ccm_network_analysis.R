library(tidyverse)

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