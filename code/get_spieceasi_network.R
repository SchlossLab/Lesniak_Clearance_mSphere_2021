
#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
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

# create output dir if doesnt already exist
for(i in c('data/process/', 'results/figures/')){
	save_dir <- paste0(i, 'spieceasi/')
	ifelse(!dir.exists(save_dir), 
		dir.create(save_dir), 
		print(paste0(save_dir, ' directory ready')))
}

antibiotic_run <- list('amp', 'clinda', 'cef', 'metro', 'strep', 'vanc',
	c("clinda", "cef", "strep"), c("cef", "metro", "strep"), c("clinda", "cef", "metro", "strep"),
	c("clinda", "vanc", "amp", "cef", "metro", "strep"))[[run_number]]
print(paste('Running SpiecEasi network analysis on', 
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

get_se_network <- function(antibiotic){
	se_cdiff_df <- meta_file %>% 
		inner_join(clearance, by = 'mouse_id') %>% 
		filter(cdiff == T, clearance != 'colonized', day > 0, abx %in% antibiotic) %>% 
		inner_join(rel_abund, by = 'group') %>% 
		mutate(Otu00Cdiff = log10(CFU + 0.01)) %>% 
		select(contains('Otu')) %>% 
		as.matrix
	se_cdiff_df <- se_cdiff_df[ ,apply(se_cdiff_df, 2, function(x) sum(x > 0)) > 0.1*nrow(se_cdiff_df)]
	# the main SPIEC-EASI pipeline: Data transformation, sparse inverse covariance estimation and model selection
	se.mb.cdiff <- spiec.easi(se_cdiff_df, method='mb', lambda.min.ratio=1e-2,
		nlambda=20, pulsar.params=list(rep.num=999))
	#for actual run, use rep.num (99 or 999) 
	## set size of vertex proportional to clr-mean
	vsize    <- rowMeans(clr(se_cdiff_df, 1))+6

	ig.mb <- getRefit(se.mb.cdiff)
	colnames(ig.mb) <- rownames(ig.mb) <- gsub('Otu0*', '', colnames(se_cdiff_df))
	# save output from spiec-easi, but a sparse matrix so need to convert to dataframe before writing
	write.table(data.frame(otu = colnames(ig.mb)), paste0('data/process/spieceasi/se_df_', 
			paste(antibiotic, collapse = '_'), '_features.txt'), 
		sep = '\t', quote = F, row.names = T)
	write.table(Matrix::summary(ig.mb), paste0('data/process/spieceasi/se_df_', 
			paste(antibiotic, collapse = '_'), '.txt'), 
		sep = '\t', quote = F, row.names = F)
		## to read in sparse matrix
		# dd <- read.table(file_name, header=TRUE)
		## has columns (i, j, x) -> we can use via do.call() as arguments to sparseMatrix():
		# mm <- do.call(sparseMatrix, dd)
	return(ig.mb)
	#cdiff.coord <- layout.fruchterman.reingold(ig.mb)
	#network_plot <- plot(ig.mb, layout=cdiff.coord, vertex.size=vsize, main = antibiotic)
	#ggsave(paste0('cdiff_network_', antibiotic, '_plot.jpg'), network_plot, width = 10, height = 10)
}

se_network_matrix <- get_se_network(antibiotic_run)

network_labels <- tax_df %>% 
	mutate(otu_number = gsub('OTU ', '', otu_label),
		tax_otu_label = gsub('_unclassified', '', tax_otu_label),
		tax_otu_label = gsub(' \\(', '\\\n\\(', tax_otu_label)) %>% 
	pull(tax_otu_label)

get_first_order <- function(data_input){
	variable_name <- antibiotic_run
	#variable_name <- gsub('_se', '', deparse(substitute(data_input)))
	print(variable_name)
	if(sum(data_input[,'Cdiff']) > 0){
		first_order_otus <- c(names(which(data_input[,'Cdiff'] > 0)), 'Cdiff')
		#second_order_otus <- names(apply(data_input[,otus], 1 , sum) > 0)
		data_input <- data_input[first_order_otus, first_order_otus]
		labels <- network_labels[as.numeric(colnames(data_input))]
		labels[length(labels)] <- 'C. difficile'
		colnames(data_input) <- rownames(data_input) <- labels
		#ig.mb <- ig.mb[second_order_otus, second_order_otus]
		first_order_network <- adj2igraph(data_input, 
			vertex.attr = list(name = colnames(data_input)))
		pdf(paste0('results/figures/spieceasi/se_network_', variable_name, '.pdf'))
			plot(first_order_network, main = paste('Spiec-easi Network for C. difficile in',
				paste(variable_name, collapse = '_')) )
		dev.off()
		} else {
			print(paste0('No C_difficile interactions for ', paste(antibiotic_run, collapse = ' ')))
	}

}

get_first_order(se_network_matrix)
