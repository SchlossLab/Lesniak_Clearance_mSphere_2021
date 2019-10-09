# using http://www.bnlearn.com/examples/useR19-tutorial/
# install.packages('bnlearn')
library(bnlearn)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("Rgraphviz")
library(Rgraphviz)
library(tidyverse)


run_number <- as.numeric(commandArgs(TRUE))

# file names relative to code directory for Rmd
meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
sum_taxa_function <- 'code/sum_otu_by_taxa.R'

# read in data
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	mutate(mouse_id = paste(cage, mouse, sep = '_'))
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F, row.names = "Group") %>% 
	select(-label, -numOtus)
tax_df <- read.table(tax_file, sep = '\t', header = T, stringsAsFactors = F)
source(sum_taxa_function) # function to create taxanomic labels for OTUs

antibiotic_run <- list('amp', 'clinda', 'cef', 'metro', 'strep', 'vanc',
		c("clinda", "vanc", "amp", "cef", "metro", "strep"))[[run_number]]
print(paste('Running bayesian network analysis on', 
	paste(antibiotic_run, collapse = ' ')))
# get relative abundances
n_seqs <- unique(apply(shared_file, 1, sum))
rel_abund <- data.frame(group = row.names(shared_file), 100*shared_file/n_seqs,
	stringsAsFactors = F)
min_rel_abund <- 100 * 1/n_seqs

# create a dataframe with the change in C. difficile CFU by day
differenced_cfu <- meta_file %>% 
	group_by(mouse_id) %>% 
	arrange(mouse_id, day) %>% 
	mutate(delta_cfu = CFU - lag(CFU)) %>% 
	ungroup %>% 
	mutate(delta_trend = ifelse(sign(delta_cfu) == -1, -1, 1),
		diff_cfu = delta_trend * log10(abs(delta_cfu)),
		diff_cfu = ifelse(diff_cfu == -Inf, 0, diff_cfu))


# create a dataframe with the change in relative abundance by day
diff_rel_abund <- meta_file %>% 
	select(group, mouse_id, day) %>% 
	inner_join(rel_abund, by = 'group') %>% 
	group_by(mouse_id) %>% 
	arrange(mouse_id, day) %>% 
	mutate_at(vars(contains('Otu')), 
		function(otu) { otu - lag(otu) } ) %>% 
	ungroup

# categorize mice based on the colonization levels at the end of the experiment
clearance <- differenced_cfu %>% 
	group_by(mouse_id) %>% 
	summarise(last_sample = max(day),
		max_cfu = max(CFU)) %>% 
	left_join(differenced_cfu, by = c('mouse_id')) %>% 
	filter(day == last_sample) %>% 
	select(mouse_id, end_point_cfu = CFU, delta_trend, diff_cfu, last_sample, max_cfu) %>% 
	mutate(clearance = case_when(max_cfu < 1 ~ 'uncolonized',
		end_point_cfu == 0 ~ 'cleared', 
		delta_trend < 0 & end_point_cfu < 100000 ~ 'clearing', # 10^5 separates the mice 
		T ~ 'colonized'))

# create labels for the 
network_labels <- tax_df %>% 
	mutate(otu_number = gsub('OTU ', '', otu_label),
		tax_otu_label = gsub('_unclassified', '', tax_otu_label),
		tax_otu_label = gsub(' \\(', '\\\n\\(', tax_otu_label)) %>% 
	select(OTU, tax_otu_label) %>% 
	rbind(data.frame(OTU = 'C_difficile', tax_otu_label = 'C. difficile',
		stringsAsFactors = F))

get_bayesian_network <- function(antibiotic){
	bn_df <- diff_rel_abund %>% 
		filter(!is.na(Otu000001)) %>% 
		inner_join(
			select(filter(differenced_cfu, abx %in% antibiotic), 
				group, C_difficile = diff_cfu) %>% 
				filter(),
			by = c('group')) %>% 
		select(-group, -mouse_id, -day)
	bn_features <- bn_df %>% 
		gather(otu, diff_abundance) %>% 
		group_by(otu) %>% 
		summarise(sd = sd(diff_abundance)) %>% 
		filter(sd > min_rel_abund)  %>% 
		pull(otu)
	bn_df <- bn_df %>% 
		select(one_of(bn_features))

	# learn bayesian network structure of directed acyclic graph
	##bn_dag <- hc(bn_df)
	## visualize dag
	##graphviz.plot(bn_dag, shape = "ellipse")
	# since data not likely all normal or linear, bootstrap to find 
	str_bn <- boot.strength(bn_df, R = 200, algorithm = "hc"#,
		#algorithm.args = list(whitelist = wl, blacklist = bl)
		)
	#plot(str_bn)
	#head(str_bn) # strength and direction of each interaction
	#attr(str_bn, "threshold") # threshold used to determine included vertices
	#avg_bn <- averaged.network(str_bn)
	#graphviz.plot(avg_bn)
	cdiff_bn <- str_bn %>% 
		.[(.$from == 'C_difficile' | .$to == 'C_difficile') & .$strength > 0.5, ]
	tax_labels <- cdiff_bn %>% 
		left_join(network_labels, by = c('from' = 'OTU')) %>% 
		left_join(network_labels, by = c('to' = 'OTU')) %>% 
		mutate(from = tax_otu_label.x,
			to = tax_otu_label.y)
	cdiff_bn$from <- tax_labels$from
	cdiff_bn$to <- tax_labels$to
	avg_bn <- averaged.network(cdiff_bn)

	pdf(paste0('results/figures/', antibiotic, '_bayesian_network_plot.pdf'))
	strength.plot(avg_bn, cdiff_bn, shape = "ellipse")
	dev.off()

	return(data.frame(str_bn, antibiotic = paste(antibiotic, collapse = '_'), stringsAsFactors = F))
}
print('Beginning bayesian network analysis')
bn_df <- get_bayesian_network(antibiotic_run)

write.table(bn_df, paste0('data/process/bn_df_', 
		paste(antibiotic_run, collapse = '_'), '.txt'), 
	sep = '\t', quote = F, row.names = F)


