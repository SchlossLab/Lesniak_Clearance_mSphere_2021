library(tidyverse)
library(igraph)
library(patchwork)
library(grid)
library(rsvg)

save_dir <- 'scratch/ccm_networks_by_genus/'
data_path <- "scratch/ccm/"   # path to the data
treatment_list <- list.files(data_path)

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

# to do
#	plot only mice used for ccm (ones with all days)
#	edit network to show strength

print(paste0('Running analysis on output from ', paste(treatment_list, collapse = ', ')))

lapply(treatment_list, function(current_treatment){
	# get file name
	# current_treatment <- 'clinda_10_FALSE'
	print(paste('Running', current_treatment))
	file <- file.path(current_treatment, 
		dir(paste0(data_path, current_treatment), pattern = "ccm_by*"))
	
	if(length(file) == 0){
		paste('No output files for ', current_treatment)
	} else {

		input_data <- read.table(file.path(data_path, file), header = T, stringsAsFactors = F)
		
		interaction_data <- input_data %>% 
			group_by(driver_otu, driven_otu) %>% 
			summarise(p_value = mean(ccm_p_value),
				strength = mean(driver_predicts_driven)) %>% 
			ungroup() %>% 
			filter(strength > mean(strength)) %>% 
			mutate(weight = case_when(p_value > 0.05 ~ 0,
				driver_otu == driven_otu ~ 0,
				TRUE ~ abs(strength))) %>% 
		#	left_join(select(taxonomic_labels, otu, tax_otu_label)) %>% 
		#	left_join(select(taxonomic_labels, otu, tax_otu_label), by = c('affected_otu'='otu')) %>% 
			select(driver_taxa = driver_otu, #tax_otu_label.x, 
				driven_taxa = driven_otu, #tax_otu_label.y, 
				p_value, strength, weight)

		interaction_matrix_plot <- interaction_data %>% 
			ggplot(aes(driver_taxa, driven_taxa)) + geom_tile(aes(fill = weight)) + 
				scale_fill_gradient(low = 'white', high = 'blue') + 
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) + 
				labs(title = paste0('Interactions from ', 
					current_treatment),
					x = 'Driver OTU', y = 'Driven OTU')

		if(sum(interaction_data$weight) > 0){
			nodes <- interaction_data %>% 
				distinct(driver_taxa) %>% 
				rename(label = driver_taxa) %>% 
				full_join(interaction_data %>% 
					distinct(driven_taxa) %>% 
					rename(label = driven_taxa),
					by = 'label') %>% 
				rowid_to_column('id')
			edges <- interaction_data %>% 
				left_join(nodes, by = c('driver_taxa' = 'label')) %>% 
				rename(from = id) %>% 
				left_join(nodes, by = c('driven_taxa' = 'label')) %>% 
				rename(to = id) %>% 
				filter(weight > 0) %>% 
				select(from, to, weight)
			# igraph		
			routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
			#network_plot <- 
			svg(paste0('scratch/temp/', current_treatment, '_temp.svg'))
			plot(routes_igraph, layout = layout_with_graphopt, 
				edge.width = E(routes_igraph)$weight*5, edge.curved = T,
				edge.arrow.size = .6)
			dev.off()
			
			img <- rsvg(paste0('scratch/temp/', current_treatment, '_temp.svg'), width = 2000)
			g <- rasterGrob(img, interpolate=TRUE, width = 1.1)
			network_plot <- ggplot() +
				annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
				geom_point()

			interacting_otus <- nodes$label

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
						facet_wrap(~bacteria, scales = 'free_y', ncol = 2) +
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
				(otu_temporal_plot | cdiff_temporal_plot) / 
				(network_plot | interaction_matrix_plot ),
				width = 20, height = 20)
		} else {
			ggsave(paste0(save_dir, current_treatment, '_dynamics_nonsig.jpg'),
				interaction_matrix_plot,
				width = 20, height = 20)
		}
	}
})