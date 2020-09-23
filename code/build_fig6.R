##############
#
# run script to generate plots for Figure 5
#	What interactions associate with clearance of C. difficile colonization?
# 
# Nick Lesniak 04-13-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/sample.final.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#
##############

library(SpiecEasi)
library(igraph)
library(tidyverse)
library(cowplot)
#library(GGally)
# inorder to have mixed text formatting in network labels
# changed geom_text() to geom_richtext() - lines 904 and 1068 
# used trace(ggnet2, edit='nano') to insert changes into ggnet2 function
# also edited R function saved in code/R/functions/ggnet2.R
library(network)
library(sna)
library(intergraph)
library(ggtext)
library(scales)
source('code/R/functions/ggnet2.R')

seed <- 18
meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/sample.final.shared'
subsampled_shared_file <- 'data/mothur/sample.final.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))
# arguments for spiec easi, 
se_pargs <- list(rep.num=99, seed=seed, ncores=4)
# function to set offset angle to arrange node labels outside circle
#  kjhealy/polar-labels.r https://gist.github.com/kjhealy/834774
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
# read in data
tax_df <- read_tsv(tax_file)
network_labels <- tax_df %>% 
	mutate(otu_number = gsub('OTU ', '', otu_label),
		tax_otu_label = gsub('_unclassified', '', tax_otu_label), # remove unclassified note
		tax_otu_label = paste0('*', tax_otu_label), # italicize
		tax_otu_label = gsub(' \\(', '*<br/>(', tax_otu_label), # add line breaks between tax and otu
		tax_otu_label = gsub('_', ' ', tax_otu_label), #convert underscores to spaces
		tax_otu_label = gsub(' ([XI1V]+)\\*<', '\\* \\1<', tax_otu_label), # unitalicize roman numerals
		tax_otu_label = case_when(grepl('OTU 10)|OTU 45)', tax_otu_label) ~ # bold otus in multiple networks
			paste0('**', tax_otu_label, '**'),
			T ~ tax_otu_label)) %>% 
	pull(tax_otu_label)

meta_df   <- read_tsv(meta_file) %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'),
		cdiff == T, clearance != 'uncolonized', day > 0,)
subsampled_shared_df <- read_tsv(subsampled_shared_file) %>% 
	select(-numOtus, -label) %>% 
	filter(Group %in% meta_df$group)
shared_df <- read_tsv(shared_file) %>% 
	select(-numOtus, -label) %>% 
	filter(Group %in% meta_df$group)

get_cdiff_network <- function(antibiotic, clearance_status){
	# select color for antibiotic
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)

	se_df <- meta_df %>% 
		filter(abx == antibiotic,
			clearance %in% clearance_status) %>% 
		select(group, Cdiff = log10CFU) %>% 
		filter(!is.na(Cdiff)) %>% 
		inner_join(shared_df, by = c('group' = 'Group')) %>% 
		select(-group)

	# select only OTUs present in 10% of samples
	otus_present <- se_df %>% 
		summarise_all(function(x){
				sum(x >= 1) >= .1 * nrow(shared_df)
			}) %>% 
		gather(OTU, present) %>% 
		filter(present == T) %>% 
		pull(OTU)
	se_df <- se_df %>% 
		select(otus_present) %>% 
		as.matrix

	# SPIEC-EASI: data transformation, sparse inverse covariance estimation and model selection
	se_model <- spiec.easi(se_df, method = 'mb', lambda.min.ratio = 1e-3, nlambda = 500,
		sel.criterion = 'bstars', pulsar.select = TRUE, pulsar.params = se_pargs)
	
	# set size of vertex proportional to clr-mean
	vsize    <- rowMeans(clr(se_df, 1))+6
	# determine edge weights
	se_beta <- symBeta(getOptBeta(se_model), mode='maxabs')
	se_edges <- Matrix::summary(se_beta)
	# network degree distributions
	se_network <- adj2igraph(getRefit(se_model))
	se_dd <- degree.distribution(se_network)
	# determine network stability (closest to 0.05 is best, increase nlambda if not close)
	se_stability <- getStability(se_model)
	if(se_stability < 0.045){ stop(paste0('Stability low (', se_stability, 
		'), increase nlambda or decrease lambda.min.ratio arguments'))}

	## setup model output for graphing network
	se_interaction_matrix <- getRefit(se_model)
	# name matrix positions
	colnames(se_interaction_matrix) <- rownames(se_interaction_matrix) <- gsub('Otu0*', '', colnames(se_df))
	# subset network to only OTUs directly interecting with C. difficile
	first_order_otus <- c(names(which(se_interaction_matrix[,'Cdiff'] > 0)), 'Cdiff')
	#second_order_otus <- names(apply(se_data[,otus], 1 , sum) > 0)
	cdiff_interactions <- se_interaction_matrix[first_order_otus, first_order_otus]
	labels <- c(network_labels[as.numeric(head(colnames(cdiff_interactions), -1))], '*C. difficile*')
	colnames(cdiff_interactions) <- rownames(cdiff_interactions) <- labels
	#se_interaction_matrix <- se_interaction_matrix[second_order_otus, second_order_otus]
	# add edge weights
	# names matrix positions
	colnames(se_beta) <- rownames(se_beta) <- gsub('Otu0*', '', colnames(se_df))
	# subset network to OTUs interacting with C difficile
	wt_cdiff_interactions <- se_beta[first_order_otus, first_order_otus]
	colnames(wt_cdiff_interactions) <- rownames(wt_cdiff_interactions) <- labels
	# create igraph object with edge weights
	wt_first_order_network <- adj2igraph(wt_cdiff_interactions, 
		vertex.attr = list(name = colnames(wt_cdiff_interactions)))

	# setup network attributes to create igraph network graph for output
	vsize <- vsize[gsub('Otu0*', '', names(vsize)) %in% first_order_otus]
	edge_wt <- abs(E(wt_first_order_network)$weight)
	edge_direction <- ifelse(E(wt_first_order_network)$weight < 0, 'red', 'blue')
	lab.locs <- radian.rescale(x=1:length(first_order_otus), direction=-1, start=0)
	network <- adj2igraph(cdiff_interactions, 
		vertex.attr = list(name = colnames(cdiff_interactions), size = vsize^2/10, 
			color = abx_col, label.color='black', label.cex = 0.7, label.dist = 2, 
			label.degree = lab.locs),
		edge.attr = list(width = edge_wt*10, color = edge_direction))

	se_output <- list(edge_wts = se_edges, degree_dist = se_dd, stability = se_stability, 
		all_otus = colnames(se_df), otus = first_order_otus, full_network = se_network, 
		cdiff_network = network, antibiotic = antibiotic, 
		clearance = paste(clearance_status, collapse = '_'))
	return(se_output)
}

clinda_network <- get_cdiff_network('Clindamycin', 'Cleared')
strep_network <- get_cdiff_network('Streptomycin', 'Cleared')
cef_network <- get_cdiff_network('Cefoperazone', 'Cleared')
strep_colonized_network <- get_cdiff_network('Streptomycin', 'Colonized')
cef_colonized_network <- get_cdiff_network('Cefoperazone', 'Colonized')

set.seed(2)
clinda_network_graph <- ggnet2(clinda_network$cdiff_network, mode = 'kamadakawai',
		color = 'color', label = T, size = 'size', vjust = 1.3, label.size = 3.5,
		edge.size = 'width', edge.color = 'color', layout.exp = 0.2) +
	guides(size = FALSE) +
	ylim(-0.1, 1) + xlim(-0.1, 1.1) + 
	theme(axis.title = element_blank(), axis.text =element_blank(),
	    axis.ticks = element_blank())
cef_network_graph <- ggnet2(cef_network$cdiff_network, mode = 'kamadakawai',
		color = 'color', label = T, size = 'size', vjust = 1.3, label.size = 3.5,
		edge.size = 'width', edge.color = 'color', layout.exp = 0.2) +
	guides(size = FALSE) +
	ylim(-0.4, 1.2) + xlim(-0.25, 1.25) + 
	theme(axis.title = element_blank(), axis.text = element_blank(),
	    axis.ticks = element_blank())
strep_network_graph <- ggnet2(strep_network$cdiff_network, mode = 'kamadakawai',
		color = 'color', label = T, size = 'size', vjust = 1.3, label.size = 3.5,
		edge.size = 'width', edge.color = 'color', layout.exp = 0.2) +
	guides(size = FALSE) + 
	ylim(-0.1, 1) + xlim(-0.1, 1.1) + 
	theme(axis.title = element_blank(), axis.text = element_blank(),
	    axis.ticks = element_blank())


networks <- list(clinda_network, strep_network, cef_network, strep_colonized_network, cef_colonized_network)
# centrality
# all look fairly similar
# slightly lower amount of high degree for comminities remaining colonized
# Cef has significantly different betweenness, 
#  cleared communities have much higher betweenness centrality
#   so cleared communities have slightly more connections 
get_centrality <- function(x){
	net <- x$full_network
	tibble(antibiotic = x$antibiotic,
		clearance = x$clearance,
		#otu = x$all_otus,
		degree = igraph::degree(net, mode="in"), # number of its adjacent edges
		betweenness = igraph::betweenness(net, directed=T, weights=NA)) %>% # the number of shortest paths going through node
		gather(metric, value, -antibiotic, -clearance)
}

centrality_df <- map_dfr(networks, get_centrality) 
# test diffs
pvalue_df <- centrality_df %>% 
	mutate(subset = paste(antibiotic, clearance, sep = '_')) %>% 
	select(metric, subset, value) %>% 
	group_by(metric) %>% 
	nest()  %>% 
	mutate(data = map(data, function(data) nest(group_by(data, subset)))) %>% 	
	mutate(data = map(data, function(nested_df){
		test_df <- sapply(nested_df$data, function(x) sapply(nested_df$data, function(y) wilcox.test(x$value,y$value)$p.value)) %>% 
			data.frame
		test_df[upper.tri(test_df)] <- NA # set all of upper triangle to 1 to eliminate duplicate comparisons
		colnames(test_df) <- nested_df$subset
		test_df <- test_df %>% 
			mutate(row_names = nested_df$subset) %>% 
			pivot_longer(col = -row_names, names_to = 'col_names', values_to = 'pvalue') %>% 
			filter(pvalue != 1, !is.na(pvalue)) %>% # eliminate all self comparisons and upper triangle
			mutate(pvalue = p.adjust(pvalue, method = 'BH')) # correct p values
		return(test_df)
		})) %>% 
	unnest(data)

annotation_df <- data.frame(metric = rep(c('Degree', 'Betweenness'), each = 6),
	x1=c(1, 1, 2.85, 2.85, 3.15, 3.15, 1, 1, 2.85, 2.85, 3.15, 3.15), 
	x2=c(1.85, 2.15, 1.85, 2.15, 1.85, 2.15, 1.85, 2.15, 1.85, 2.15, 1.85, 2.15), 
	xnote = c(1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5),
	y1 = c(10^0.843, 10^0.878, 10^0.913, 10^0.948, 10^0.983, 10^1.018, 
			10^2.41, 10^2.51, 10^2.61, 10^2.71, 10^2.81, 10^2.91),
	ynote = c(10^0.8445, 10^0.8795, 10^0.9145, 10^0.9495, 10^0.9845, 10^1.0195, 
			10^2.42, 10^2.52, 10^2.62, 10^2.72, 10^2.82, 10^2.92),
	annotations='*', clearance = 'NA', antibiotic = 'NA')

centrality_plot <- centrality_df %>% 
	mutate(antibiotic = factor(antibiotic, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin')),
		metric = case_when(metric == 'betweenness' ~ 'Betweenness',
			metric == 'degree' ~ 'Degree',
			T ~ metric)) %>% 
	ggplot(aes(x = antibiotic, y = value, fill = clearance, color = antibiotic)) + 
		geom_boxplot(width = 0.6,  position = position_dodge2(preserve = "single"),
			show.legend = F) + 
		geom_point(size = NA, shape = 22) + 
		scale_color_manual(values = abx_color$color, limits = abx_color$abx) + 
		scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) + 
		facet_wrap(.~metric, scales = 'free') + 
		scale_y_log10() + 
		guides(color = 'none') + theme_bw() + 
		labs(x = NULL, y = NULL, fill = NULL) +
		theme(legend.position = c(0.5, 0.565),
			panel.spacing = unit(5,'lines'),
			strip.background = element_blank(),
			strip.text = element_text(size = 14)) +
		guides(fill = guide_legend(override.aes = list(size = 5))) + 
		geom_text(data = annotation_df, aes(x = xnote, y = ynote, label = annotations), 
			color = 'black') +
		geom_segment(data = annotation_df, aes(x = x1, xend = x2, y = y1, yend = y1), 
			color = 'black', size = 0.25)

write_tsv(pvalue_df, 'data/process/network_wilcoxon_BHadjusted_pvalue.tsv')

ggsave('results/figures/figure_6.jpg',
		plot_grid(
			plot_grid(
				plot_grid(clinda_network_graph, labels = c('Clindamycin'), label_colour = '#A40019'),
				plot_grid(plot_grid(NULL, cef_network_graph, ncol = 1, rel_heights = c(1, 10)), 
					labels = c('Cefoperazone'), label_colour = '#3A9CBC'),
				plot_grid(strep_network_graph, labels = c('Streptomycin'), label_colour = '#D37A1F'), nrow = 1),
		centrality_plot, 
		ncol = 1, rel_heights = c(3,2), labels = c('A', 'B')), 
	width = 18, height = 13)
