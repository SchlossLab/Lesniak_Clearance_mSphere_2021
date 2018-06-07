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

	otu_names <- gsub('Otu0+', '', unique(data$otu))
	otu_names[otu_names == 'CFU'] <- 'C. difficile'

	interaction_data <- data %>% 
		group_by(otu, affected_otu) %>% 
		summarise(p_value = median(p_value),
			strength = median(strength)) %>% 
		mutate(adj_strength = ifelse(p_value > 0.05, 0, strength)) 

	interaction_data %>% 
		ggplot(aes(otu, affected_otu)) + geom_tile(aes(fill = adj_strength)) + 
			scale_fill_gradient(low = 'white', high = 'blue')

test_data <- interaction_data %>% 
	filter(otu1 %in% c('Otu000001', 'Otu000003', 'Otu000004'), 
		otu2 %in% c('Otu000001', 'Otu000003', 'Otu000004'))

	num_nodes <- length(unique(data$otu))
	network_matrix <- matrix(interaction_data$adj_strength,
		nrow = num_nodes, ncol = num_nodes)
	diag(network_matrix) <- 0

	network <- as.network(x = network_matrix,
		directed = T,
		loops = F,
		matrix.type = 'adjacency')

	network.vertex.names(network) <- otu_names
	
	plot.network(network, # our network object
             vertex.col = 'gray', # color nodes by gender
             vertex.cex = 2, # size nodes by their age
             displaylabels = T, # show the node names
             label.pos = 5 # display the names directly over nodes
             )

	

}
