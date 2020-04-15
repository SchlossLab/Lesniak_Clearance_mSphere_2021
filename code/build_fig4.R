##############
#
# run script to generate plots for Figure 4
#	What interactions associate with clearance of C. difficile colonization?
# 
# Nick Lesniak 04-13-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#
##############

library(SpiecEasi)
library(igraph)
library(tidyverse)
library(cowplot)

seed <- 18
meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))
# arguments for spiec easi, 
se_pargs <- list(rep.num=99, seed=seed, ncores=4)
# read in data
tax_df <- read_tsv(tax_file)
network_labels <- tax_df %>% 
	mutate(otu_number = gsub('OTU ', '', otu_label),
		tax_otu_label = gsub('_unclassified', '', tax_otu_label),
		tax_otu_label = gsub(' \\(', '\\\n\\(', tax_otu_label)) %>% 
	pull(tax_otu_label)

meta_df   <- read_tsv(meta_file) %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'),
		cdiff == T, clearance != 'uncolonized', day > 0,)
shared_df <- read_tsv(shared_file) %>% 
	select(-numOtus, -label, -X6631) %>% 
	filter(Group %in% meta_df$group)
# select only OTUs present in 10% of samples
otus_present <- shared_df %>% 
	select(-Group) %>% 
	summarise_all(function(x){
			sum(x >= 1) >= .1 * nrow(shared_df)
		}) %>% 
	gather(OTU, present) %>% 
	filter(present == T) %>% 
	pull(OTU)
shared_df <- shared_df %>% 
	select(Group, one_of(otus_present)) %>% 
	left_join(select(meta_df, group, Cdiff = log10CFU), by = c('Group' = 'group'))


get_cdiff_network <- function(antibiotic, clearance_status){
	# select color for antibiotic
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)

	se_df <- meta_df %>% 
		filter(abx == antibiotic,
			clearance %in% clearance_status) %>% 
		select(group) %>% 
		left_join(shared_df, by = c('group' = 'Group')) %>% 
		filter(!is.na(Cdiff)) %>% 
		select(-group) %>% 
		as.matrix
	# SPIEC-EASI: data transformation, sparse inverse covariance estimation and model selection
	se_model <- spiec.easi(se_df, method = 'mb', lambda.min.ratio = 1e-3, nlambda = 100,
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
	labels <- c(network_labels[as.numeric(head(colnames(cdiff_interactions), -1))], 'C. difficile')
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
			color = abx_col, label.color='black', label.size = 8, label.dist = 2),
		edge.attr = list(width = edge_wt*10, color = edge_direction))

	se_output <- list(edge_wts = se_edges, degree_dist = se_dd, stability = se_stability, 
		all_otus = colnames(se_df), otus = first_order_otus, full_network = se_network, 
		cdiff_network = network)
	return(se_output)
}

clinda_network <- get_cdiff_network('Clindamycin', 'Cleared')
strep_network <- get_cdiff_network('Streptomycin')
cef_network <- get_cdiff_network('Cefoperazone')

# test for differences in centrality
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5927471/
# https://kateto.net/netscix2016.html
#bind_rows(tibble(abx = 'Clindamycin',
#		local_centr = centr_degree(clinda_network$full_network)$res, 
#		closeness = closeness(clinda_network$full_network),
#		otu = clinda_network$all_otus),
#	tibble(abx = 'Cefoperazone',
#		local_centr = centr_degree(cef_network$full_network)$res, 
#		closeness = closeness(cef_network$full_network),
#		otu = cef_network$all_otus),
#	tibble(abx = 'Streptomycin',
#		local_centr = centr_degree(strep_network$full_network)$res, 
#		closeness = closeness(strep_network$full_network),
#		otu = strep_network$all_otus))
#
# no significant difference in degree distribution
#bind_rows(cbind(abx = 'Clindamycin', as.data.frame(clinda_network$edge_wts), stringsAsFactors = F),
#	cbind(abx = 'Cefoperazone', as.data.frame(cef_network$edge_wts), stringsAsFactors = F),
#	cbind(abx = 'Streptomycin', as.data.frame(strep_network$edge_wts), stringsAsFactors = F)) %>% 
#	ggplot(aes(x, fill = abx)) + 
#		geom_histogram(binwidth = 0.1) + 
#		facet_wrap(abx~., scales = 'free_y')
# no significant difference in degree distribution
#bind_rows(tibble(abx = 'Clindamycin', Frequency = clinda_network$degree_dist, 
#		Degree = 0:(length(clinda_network$degree_dist) - 1)),
#	tibble(abx = 'Cefoperazone', Frequency = cef_network$degree_dist, 
#		Degree = 0:(length(cef_network$degree_dist) - 1)),
#	tibble(abx = 'Streptomycin', Frequency = strep_network$degree_dist, 
#		Degree = 0:(length(strep_network$degree_dist) - 1))) %>% 
#	ggplot(aes(x = Degree, y = Frequency, color = abx)) + 
#	geom_line()


	jpeg(paste0('../results/figures/spieceasi/se_', 
		antibiotic, '_network.jpg'),
		width = 5, height = 5, units = 'in', res = 150)
#	dev.off()
}

