library(tidyverse)
library(statnet)
library(geomnet)
library(patchwork)

save_dir <- 'scratch/ccm_networks_by_genus/'
data_path <- "scratch/ccm_all/"   # path to the data
files <- dir(data_path, pattern = "ccm_by*") # get file names
treatment_list <- unique(gsub('\\d{,2}.txt', '', files))
source('code/taxa_labels.R')
taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_df   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) 
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_df <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F)
source('code/sum_otu_by_taxa.R')
taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'
shared_by_genus <- sum_otu_by_taxa(taxonomy_file = taxonomy_file, 
	otu_df = shared_df, 
	taxa_level = 'genus')

for(treatment in treatment_list){
	current_treatment <- gsub('.*_data_(.+)_seed', '\\1', treatment)
	files <- dir(data_path, pattern = paste0(treatment, "*")) # get file names
	input_data <- files %>%
	  # read in all the files, appending the path before the filename
	  map(~ read.table(file.path(data_path, .), header = T, stringsAsFactors = F)) %>% 
	  reduce(rbind)

	input_data <- bind_rows(
		select(input_data, otu = otu1, strength = otu1_cause_otu2, p_value = pval_otu1_casue_otu2, affected_otu = otu2,
			prediction_slope = otu1_prediction_slope, p_slope = otu1_prediction_slope_p, E = E_A),
		select(input_data, otu = otu2, strength = otu2_cause_otu1, p_value = pval_otu2_casue_otu1, affected_otu = otu1,
			prediction_slope = otu2_prediction_slope, p_slope = otu2_prediction_slope_p, E = E_B))

#	taxonomic_labels <- get_taxa_labels(taxa_file = taxonomy_file,  taxa_level='genus', 
#		otu_subset=unique(input_data$otu)) %>% 
#		bind_rows(data.frame(otu = 'CFU', taxa = 'C. difficile', otu_label = 'C. difficile', 
#			tax_otu_label = 'C. difficile', stringsAsFactors = F))
	#otu_names <- unique(input_data$otu)
	#otu_names[otu_names == 'CFU'] <- 'C. difficile'

	interaction_data <- input_data %>% 
		group_by(otu, affected_otu) %>% 
		summarise(p_value = median(p_value),
			strength = median(strength)) %>% 
		ungroup() %>% 
		mutate(adj_strength = ifelse(p_value > 0.05, 0, strength),
			adj_strength = ifelse(otu == affected_otu, 0, adj_strength)) %>% 
#		left_join(select(taxonomic_labels, otu, tax_otu_label)) %>% 
#		left_join(select(taxonomic_labels, otu, tax_otu_label), by = c('affected_otu'='otu')) %>% 
		select(driver_taxa = otu, #tax_otu_label.x, 
			driven_taxa = affected_otu, #tax_otu_label.y, 
			p_value, strength, adj_strength)

	ggsave(paste0(save_dir, current_treatment, '_interaction_matrix.jpg'),
		interaction_data %>% 
			ggplot(aes(driver_taxa, driven_taxa)) + geom_tile(aes(fill = adj_strength)) + 
				scale_fill_gradient(low = 'white', high = 'blue') + 
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) + 
				labs(title = paste0('Interactions from ', 
					gsub('([[:alpha:]]+_){3}|_[[:alpha:]]+$', '', treatment)),
					x = 'Driver OTU', y = 'Driven OTU'),
		width = 12, height = 12)
	if(sum(interaction_data$adj_strength) > 0){
		num_nodes <- length(unique(input_data$otu))
		network_matrix <- matrix(interaction_data$adj_strength,
			nrow = num_nodes, ncol = num_nodes)
		diag(network_matrix) <- 0

		network <- as.network(x = network_matrix,
			directed = T,
			loops = F,
			matrix.type = 'adjacency')

		network.vertex.names(network) <- interaction_data$driver_taxa

		#plot.network(network, # our network object
	    #         vertex.col = 'gray', # color nodes by gender
	    #         vertex.cex = 2, # size nodes by their age
	    #         displaylabels = T, # show the node names
	    #         label.pos = 5 # display the names directly over nodes
	    #         )

	    interacting_otus <- interaction_data %>% 
	    	filter(adj_strength > 0) %>% 
	    	gather(role, taxa, driver_taxa, driven_taxa) %>% 
	    	pull(taxa) %>% 
	    	unique
		gg_network_data <- data.frame(from_id = interaction_data$driver_taxa,
			to_id = interaction_data$driven_taxa) %>% 
			mutate(to_id = ifelse(interaction_data$adj_strength == 0, NA, interaction_data$driven_taxa)) %>% 
			filter(from_id %in% interacting_otus) %>% 
			mutate(from_id = gsub(' \\(', '\n(', from_id),
				to_id = gsub(' \\(', '\n(', to_id),)

		## Using Name1 as the from node column and Name2 as the to node column.
		## If this is not correct, rewrite dat so that the first 2 columns are from and to node, respectively.

		## Joining edge and node information by from_id and label respectively.

		# create plot
		set.seed(1)
		ggsave(paste0(save_dir, current_treatment, '_network.jpg'),
			ggplot(data = gg_network_data, aes(from_id = from_id, to_id = to_id)) +
				geom_net(layout.alg = "kamadakawai", 
					size = 2, labelon = TRUE, vjust = -0.6, ecolour = "grey60",
					directed = TRUE, fontsize = 3, ealpha = 0.5) +
				xlim(c(-0.05, 1.05)) +
				theme_net() +
				theme(legend.position = "bottom"),
			width = 10, height = 10)


		#interacting_otus <- taxonomic_labels %>% 
		#	filter(tax_otu_label %in% interacting_otus) %>% 
		#	pull(otu) %>% 
		#	c(., 'CFU')
		interacting_otus <- unique(gg_network_data$from_id)

		abx_df <- meta_df %>% 
			select(group, cage, mouse, day, C_difficile=CFU, cdiff, abx, dose, delayed) %>% 
			unite(treatment, abx, dose, delayed) %>% 
			filter(cdiff == T, day >= 0, treatment == current_treatment) %>% 
			mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
			inner_join(select(shared_by_genus, one_of(c(interacting_otus, 'Group'))),
				by = c('group' = "Group"))


		otu_temporal_plot <- abx_df %>% 
			select(cage, mouse, day, one_of(interacting_otus)) %>% 
			gather(bacteria, counts, one_of(interacting_otus)) %>% 
				ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), 
					group = interaction(cage, mouse))) + 
					geom_line() + 
					facet_grid(bacteria~., scales = 'free_y') +
					theme_bw() + 
					labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
						title = 'Temporal Dynamics', subtitle = 'Colored by mouse') + 
					scale_x_continuous(breaks=seq(0,10, 1)) + 
					theme_bw(base_size = 8) + 
					theme(legend.position = 'none')
		cdiff_temporal_plot <- abx_df %>% 
			select(cage, mouse, day, C_difficile) %>% 
				ggplot(aes(x = day, y = C_difficile + 1, color = interaction(as.factor(mouse), as.factor(cage)), 
					group = interaction(cage, mouse))) + 
					geom_line() + scale_y_log10() + 
					theme_bw() + 
					labs(x = 'Day', y = 'Abundance \n (C difficle = CFU)', 
						title = 'Temporal Dynamics', subtitle = 'Colored by mouse') + 
					scale_x_continuous(breaks=seq(0,10, 1)) + 
					theme_bw(base_size = 8) + 
					theme(legend.position = 'none')


		ggsave(paste0(save_dir, current_treatment, '_dynamics.jpg'),
			otu_temporal_plot + cdiff_temporal_plot + plot_layout(ncol = 1),
			width = 10, height = 20)
	}

}
