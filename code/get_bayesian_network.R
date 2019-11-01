# using http://www.bnlearn.com/examples/useR19-tutorial/
# install.packages('bnlearn')
library(bnlearn)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("Rgraphviz")
library(Rgraphviz)
library(tidyverse)
library(forecast)


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
	save_dir <- paste0(i, 'bayesnet/')
	ifelse(!dir.exists(save_dir), 
		dir.create(save_dir), 
		print(paste0(save_dir, ' directory ready')))
}

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
		filter(day > 0) %>% 
		select(-group, -day)
	bn_features <- bn_df %>% 
		gather(otu, diff_abundance) %>% 
		group_by(otu) %>% 
		summarise(sd = sd(diff_abundance)) %>% 
		filter(sd > min_rel_abund)  %>% 
		pull(otu)
	mouse_list <- unique(bn_df$mouse_id)
	for(i in mouse_list){
		training.set <- filter(bn_df, mouse_id != i) %>% 
			select(one_of(bn_features))
		test.set <- filter(bn_df, mouse_id == i) %>% 
			select(one_of(bn_features))
		train_bn <- boot.strength(training.set, R = 10, algorithm = "hc") # learn BN structure on training set data 
		avg.diff <- averaged.network(train_bn)
		strength.plot(avg.diff, train_bn)
		fit_bn <- bn.fit(averaged.network(train_bn), training.set) # learning of parameters
		pred <- predict(fit_bn, "C_difficile", test.set) # predicts the value of node C given test set
		pred_df <- cbind(pred, test.set[, "C_difficile"], mouse_id = i) # compare the actual and predicted
		accuracy(f = pred, x = test.set$C_difficile)
	
xval <- bn.cv(training.set, bn = "hc",
	  loss = "cor-lw",
	  loss.args = list(target = "C_difficile", n = 200), runs = 10)
	
err <- numeric(10)
	
for (i in 1:10) {

	  tt <- table(unlist(sapply(xval[[i]], '[[', "observed")),
	             unlist(sapply(xval[[i]], '[[', "predicted")) > 0.50)

	  err[i] <- (sum(tt) - sum(diag(tt))) / sum(tt)

	}#FOR
	
summary(err)



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
		.[(.$from == 'C_difficile' | .$to == 'C_difficile') & 
			.$strength > attr(str_bn, "threshold"), ]
	tax_labels <- cdiff_bn %>% 
		left_join(network_labels, by = c('from' = 'OTU')) %>% 
		left_join(network_labels, by = c('to' = 'OTU')) %>% 
		mutate(from = tax_otu_label.x,
			to = tax_otu_label.y)
	cdiff_bn$from <- tax_labels$from
	cdiff_bn$to <- tax_labels$to
	avg_bn <- averaged.network(cdiff_bn)

	pdf(paste0('results/figures/bayesnet/bayesian_network_', antibiotic, '_plot.pdf'))
	strength.plot(avg_bn, cdiff_bn, shape = "ellipse")
	dev.off()

	return(data.frame(str_bn, antibiotic = paste(antibiotic, collapse = '_'), stringsAsFactors = F))
}
print('Beginning bayesian network analysis')
bn_df <- get_bayesian_network(antibiotic_run)

write.table(bn_df, paste0('data/process/bayesnet/bn_df_', 
		paste(antibiotic_run, collapse = '_'), '.txt'), 
	sep = '\t', quote = F, row.names = F)

# compare cleared and colonized networks
	bn_df <- diff_rel_abund %>% 
		filter(!is.na(Otu000001)) %>% 
		inner_join(
			select(filter(differenced_cfu, abx == antibiotic), 
				group, C_difficile = diff_cfu),
			by = c('group')) %>% 
		inner_join(select(clearance, mouse_id, clearance) %>% 
			unique, by = c('mouse_id')) %>% 
		select(-group, -mouse_id, -day) 
	bn_features <- bn_df %>% 
		select(-clearance) %>% 
		gather(otu, diff_abundance) %>% 
		group_by(otu) %>% 
		summarise(sd = sd(diff_abundance)) %>% 
		filter(sd > min_rel_abund)  %>% 
		pull(otu)
	bn_df <- bn_df %>% 
		select(clearance, one_of(bn_features))

bn_cleared <- bn_df %>% 
	filter(clearance == 'cleared') %>% 
	select(-clearance)
bn_colonized <- bn_df %>% 
	filter(clearance == 'colonized') %>% 
	select(-clearance)
# learn bayesian network structure of directed acyclic graph
# since data not likely all normal or linear, bootstrap to find 
str_colonized_bn <- boot.strength(bn_colonized, R = 200, algorithm = "hc")
str_cleared_bn <- boot.strength(bn_cleared, R = 200, algorithm = "hc")
str_bn <- boot.strength(select(bn_df, -clearance), R = 200, algorithm = "hc")

plot_bn <- function(bn_data){
	cdiff_bn <- bn_data %>% 
		.[(.$from == 'C_difficile' | .$to == 'C_difficile') & 
			.$strength > attr(bn_data, "threshold"), ]
	tax_labels <- cdiff_bn %>% 
		left_join(network_labels, by = c('from' = 'OTU')) %>% 
		left_join(network_labels, by = c('to' = 'OTU')) %>% 
		mutate(from = tax_otu_label.x,
			to = tax_otu_label.y)
	cdiff_bn$from <- tax_labels$from
	cdiff_bn$to <- tax_labels$to
	avg_bn <- averaged.network(cdiff_bn)
	strength.plot(avg_bn, cdiff_bn, shape = 'ellipse')
}
plot_bn(str_colonized_bn)
plot_bn(str_cleared_bn)
plot_bn(str_bn)

cdiff_otus <- str_bn %>% 
	filter(strength > attr(str_bn, "threshold"),
		from == 'C_difficile') %>% 
	pull(to)

cdiff_bn_subset <- str_bn %>% 
	.[.$strength > 0.5, ] %>% 
	.[.$from %in% c('C_difficile', cdiff_otus) |
		.$to %in% c('C_difficile', cdiff_otus), ]
avg_bn_subset <- averaged.network(cdiff_bn_subset)
strength.plot(avg_bn_subset, cdiff_bn_subset)

fitted <- bn.fit(avg_bn_subset, select(bn_df, one_of(unique(cdiff_bn_subset$to))))


xval <- bn.cv(bn_colonized, bn = "hc", loss = "cor-lw",
		loss.args = list(target = 'C_difficile', n = 200), runs = 10)

predcor[var] = mean(sapply(xval, function(x) attr(x, "mean")))

}#FOR

round(predcor, digits = 3)

dCoGo dGoPg dIMPA  dCoA dPPPM  dANB 
0.850 0.904 0.233 0.923 0.410 0.643


plot_diff_corr <- function(bn_data, antibiotic){
	meta_file %>% 
		inner_join(rel_abund, by = 'group') %>% 
		select(mouse_id, CFU, day, abx, 
			one_of(cdiff_bn$from)) %>% 
		filter(day >= 0, abx %in% antibiotic) %>% 
		mutate(CFU = log10(CFU + 60)) %>% 
		mutate_at(vars(contains('Otu')), function(otu) {log10(otu + 0.04)}) %>% 
		gather(otu, rel_abund, CFU, contains('Otu00')) %>% 
		left_join(select(clearance, mouse_id, clearance)) %>% 
		#left_join(tax_df, by = c('otu' = 'OTU')) %>% 
		ggplot(aes(y = rel_abund, x = day)) + 

			geom_line(aes(color = otu, group = mouse_id)) + 
			theme_bw() + 
			facet_grid(otu~clearance, scales = 'free_y') + 
			labs(x = 'Relative Abundance', y = 'C.difficile CFU')
}
diff_rel_abund %>% 
	filter(!is.na(Otu000001)) %>% 
	gather(OTU, abundance, one_of(cdiff_bn$from)) %>% 
	select(group, mouse_id, OTU, abundance) %>% 
	inner_join(
		select(filter(differenced_cfu, abx == antibiotic), 
			group, C_difficile = diff_cfu),
		by = c('group')) %>% 
	mutate(abundance = case_when(abundance > 0 ~ log10(abundance),
			abundance < 0 ~ -log10(-abundance),
			abundance == 0 ~ 0)) %>% 
	inner_join(select(clearance, mouse_id, clearance) %>% 
		unique, by = c('mouse_id')) %>% 
	filter(clearance %in% c('colonized', 'cleared')) %>% 
	ggplot( aes(x = abundance, y = C_difficile, color = OTU)) + 
		geom_point(alpha = 0.2) + 
		facet_grid(OTU~clearance, scales = 'free_x') + 
		theme_bw()