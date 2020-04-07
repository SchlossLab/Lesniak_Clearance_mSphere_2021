##############
#
# run script to generate plots for Figure 2
#	Can other conditions naturally clear C. difficile colonization?
# 
# Nick Lesniak 04-06-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#	data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.ave-std.summary
#	data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist
#	code/read.dist.R
#
##############


library(tidyverse)
library(cowplot)


meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
sum_taxa_function <- 'code/sum_otu_by_taxa.R'
alpha_div_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.ave-std.summary'
beta_div_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist'
dist_function <- 'code/read.dist.R'

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', abx_col))

# read in data
meta_df   <- read_tsv(meta_file) %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'))
shared_df <- read_tsv(shared_file) %>% 
	select(-label, -numOtus) %>% 
	filter(Group %in% meta_df$group)
tax_df <- read_tsv(tax_file)
alpha_df <- read_tsv(alpha_div_file) %>% 
	filter(group %in% meta_df$group)
source(sum_taxa_function) # function to create taxanomic labels for OTUs
	# sum_otu_by_taxa(taxonomy_df, otu_df, taxa_level = 'NA', top_n = 0, silent = T){
source(dist_function) # function to read in distance file and convert from triangle to dataframe
beta_df <- read_dist(beta_div_file) %>% 
	filter(rows %in% meta_df$group,
		columns %in% meta_df$group)

plot_diversity <- function(antibiotic){
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)
	abx_meta <- meta_df %>% 
		filter(abx == antibiotic,
			cdiff == T) %>% 
		select(group, mouse_id, day, dose) 

	alpha_sobs_plot <- alpha_df %>% 
		select(group, sobs) %>% 
		inner_join(abx_meta, by = c('group')) %>% 
		ggplot(aes(x = day, y = sobs)) + 
			geom_violin(aes(group = cut_width(day, 1)), scale = 'width', fill = abx_col, color = abx_col) + 
	        scale_x_continuous(breaks = -1:10) +
			theme_bw() + labs(x = 'Day', y = expression(~S[obs])) + 
			facet_wrap(dose~., ncol = 1)

	alpha_invsimp_plot <- alpha_df %>% 
		select(group, invsimpson) %>% 
		inner_join(abx_meta, by = c('group')) %>% 
		ggplot(aes(x = day, y = invsimpson)) + 
			geom_violin(aes(group = cut_width(day, 1)), scale = 'width', fill = abx_col, color = abx_col) + 
	        scale_x_continuous(breaks = -1:10) +
			theme_bw() + labs(x = 'Day', y = 'Inverse Simpson')+ 
			facet_wrap(dose~., ncol = 1)

	beta_plot <- beta_df %>% 
		inner_join(abx_meta, by = c('rows' = 'group')) %>% 
		inner_join(abx_meta, by = c('columns' = 'group')) %>% 
		filter(mouse_id.x == mouse_id.y, 
			dose.x == dose.y,
			day.x == min(abx_meta$day)) %>% 
		ggplot(aes(x = day.y, y = distances, group = mouse_id.x)) + 
			scale_x_continuous(breaks = -1:10) +
			coord_cartesian(ylim = c(0,1)) +
			geom_violin(aes(group = cut_width(day.y, 1)), scale = 'width', fill = abx_col, color = abx_col) + 
			theme_bw() + 
			labs(x = 'Day', y = 'Theta yc')+ 
			facet_wrap(dose.x~., ncol = 1)

	return(plot_grid(alpha_sobs_plot, alpha_invsimp_plot, beta_plot, nrow = 1))
}

ggsave('results/figures/figure_3_clindadiv.jpg', plot_diversity('Clindamycin'))
ggsave('results/figures/figure_3_strepdiv.jpg', plot_diversity('Streptomycin'))
ggsave('results/figures/figure_3_cefdiv.jpg', plot_diversity('Cefoperazone'))

plot_day0_abundance <- function(antibiotic, n_taxa){
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)
	abx_group <- meta_df %>% 
		filter(abx == antibiotic,
			day == 0) %>% 
		pull(group)
	abx_shared <- shared_df %>% 
		filter(Group %in% abx_group)

	plot <- sum_otu_by_taxa(tax_df, abx_shared, taxa_level = 'Genus', top_n = n_taxa) %>% 
		left_join(select(meta_df, group, dose), by = c('Group' = 'group')) %>% 
		group_by(Group) %>% 
		mutate(total = sum(abundance),
			relative_abundance = log10(abundance/total * 100)) %>% 
		ggplot(aes(x = Group, y =taxa, fill = relative_abundance)) + 
			geom_tile() +
			scale_fill_gradient(low="white", high=abx_col, limits = c(0,2), na.value = NA, 
				breaks = c(0, 1, 2), labels = c('', '10', '100')) + 
			theme_bw() + 
			facet_wrap(dose~., scales = 'free_x', nrow = 1) +
			labs(x = NULL, y = NULL, #title = 'Clindamycin Community',
				fill = 'Relative Abundance (%)\nColor Intesity based on Log10') + 
			theme(axis.title.x=element_blank(),
	        	axis.text.x=element_blank(),
	        	axis.ticks.x=element_blank(),
	        	axis.text.y = element_text(angle = 45),
	        	legend.position = 'bottom')
	return(plot)
}

clinda_abundance <- plot_day0_abundance('Clindamycin', 10)
cef_abundance <- plot_day0_abundance('Cefoperazone', 20)
strep_abundance <- plot_day0_abundance('Streptomycin', 12)

plot_grid(strep_cfu_plot, cef_cfu_plot)