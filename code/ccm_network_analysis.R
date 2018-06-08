library(tidyverse)
library(statnet)
library(geomnet)
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

	otu_names <- gsub('Otu0+', 'OTU_', unique(data$otu))
	otu_names[otu_names == 'CFU'] <- 'C. difficile'

	interaction_data <- data %>% 
		group_by(otu, affected_otu) %>% 
		summarise(p_value = median(p_value),
			strength = median(strength)) %>% 
		ungroup() %>% 
		mutate(adj_strength = ifelse(p_value > 0.05, 0, strength),
			adj_strength = ifelse(otu == affected_otu, 0, adj_strength),
			otu = ifelse(otu == 'CFU', 'C. difficile', gsub('Otu0+', "OTU_", otu)),
			affected_otu = ifelse(affected_otu == 'CFU', 'C. difficile', 
				gsub('Otu0+', "OTU_", affected_otu)))

	interaction_heatmap <- interaction_data %>% 
		ggplot(aes(otu, affected_otu)) + geom_tile(aes(fill = adj_strength)) + 
			scale_fill_gradient(low = 'white', high = 'blue') + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) + 
			labs(title = paste0('Interactions from ', 
				gsub('([[:alpha:]]+_){3}|_[[:alpha:]]+$', '', treatment)),
			x = 'Effector OTU', y = 'Affected OTU')

	num_nodes <- length(unique(data$otu))
	network_matrix <- matrix(interaction_data$adj_strength,
		nrow = num_nodes, ncol = num_nodes)
	diag(network_matrix) <- 0

	network <- as.network(x = network_matrix,
		directed = T,
		loops = F,
		matrix.type = 'adjacency')

	network.vertex.names(network) <- otu_names

	#plot.network(network, # our network object
    #         vertex.col = 'gray', # color nodes by gender
    #         vertex.cex = 2, # size nodes by their age
    #         displaylabels = T, # show the node names
    #         label.pos = 5 # display the names directly over nodes
    #         )

    interacting_otus <- interaction_data %>% 
    	filter(adj_strength > 0) %>% 
    	gather(role, otu, otu, affected_otu) %>% 
    	pull(otu) %>% 
    	unique
	gg_network_data <- data.frame(from_id = interaction_data$otu,
		to_id = interaction_data$affected_otu) %>% 
		mutate(to_id = ifelse(interaction_data$adj_strength == 0, NA, interaction_data$affected_otu)) %>% 
		filter(from_id %in% interacting_otus)

	## Using Name1 as the from node column and Name2 as the to node column.
	## If this is not correct, rewrite dat so that the first 2 columns are from and to node, respectively.

	## Joining edge and node information by from_id and label respectively.

	# create plot
	set.seed(1)
	network_plot <- ggplot(data = gg_network_data, aes(from_id = from_id, to_id = to_id)) +
		geom_net(layout.alg = "kamadakawai", 
			size = 2, labelon = TRUE, vjust = -0.6, ecolour = "grey60",
			directed = TRUE, fontsize = 3, ealpha = 0.5) +
		xlim(c(-0.05, 1.05)) +
		theme_net() +
		theme(legend.position = "bottom")

	ggsave(paste0('scratch/ccm_networks/ccm_network', treatment, '.jpg'),
				plot_grid(interaction_heatmap, network_plot),
		width = 14, height = 7)

}
