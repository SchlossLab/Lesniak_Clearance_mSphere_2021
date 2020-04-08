##############
#
# run script to generate plots for Figure 1
#	What occurs while C. difficile colonization is naturally cleared?
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

# read in data
meta_df   <- read_tsv(meta_file) %>% 
	filter(abx == 'Clindamycin')
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

# plot C difficile colonization level
colonization_plot <- meta_df %>% 
	mutate(CFU = case_when(CFU == 0 ~ 60, # shift 0 counts to just below limit of detection line
		T ~ CFU)) %>% 
	filter(!is.na(CFU)) %>% 
	ggplot(aes(x = day, y = CFU)) + 
		stat_summary(fun.y=median, geom="line", size = 1, color = '#A40019') + # create median line
        stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), color = '#A40019') + # with bars for IQR
        scale_x_continuous(breaks = -1:10) + # make ticks for each day
		annotate(x = -1, y = 200, geom = 'label', label = "LOD", # create a dotted line labeled LOD for limit of detection
			fill = "white", color = 'black', label.size = NA) + 
		geom_hline(yintercept = 101, linetype = 'dashed', size = 0.25) + 
		scale_y_log10(
   			breaks = scales::trans_breaks("log10", function(x) 10^x),
   			labels = scales::trans_format("log10", scales::math_format(10^.x))) + # scale y axis log10 and label 10^x
		theme_bw() + labs(x = 'Day', y = expression(italic('C. difficile')~' CFU'))

shared_genus <- sum_otu_by_taxa(tax_df, shared_df, taxa_level = 'Genus', top_n = 10) # sum at the genus level for the top 10

# plot relative abundance over time in a heatmap plot
# mice along the x axis, taxonomic classification along the y axis and color intensity by log10 relative abundance
abundance_plot <- shared_genus %>% 
	full_join(meta_df, by = c('Group' = 'group')) %>% 
	group_by(Group) %>% 
	mutate(total = sum(abundance),
		relative_abundance = log10(abundance/total * 100),
		taxa = gsub('_unclassified', '', taxa)) %>% 
	ggplot(aes(x = mouse_id, y =taxa, fill = relative_abundance)) + 
		geom_tile() +
		scale_fill_gradient(low="white", high='#A40019', limits = c(0,2), na.value = NA, 
			breaks = c(0, 1, 2), labels = c('', '10', '100')) + 
		theme_bw() + 
		labs(x = NULL, y = NULL, #title = 'Clindamycin Community',
			fill = 'Relative Abundance (%)\nColor Intesity based on Log10') + 
		theme(axis.title.x=element_blank(),
        	axis.text.x=element_blank(),
        	axis.ticks.x=element_blank(),
        	axis.text.y = element_text(angle = 45),
        	legend.position = 'bottom') +
		facet_wrap(.~day, nrow = 1)

# plot Sobs by day
alpha_sobs_plot <- alpha_df %>% 
	select(group, sobs) %>% 
	left_join(select(meta_df, group, day), by = c('group')) %>% 
	ggplot(aes(x = day, y = sobs)) + 
		geom_violin(aes(group = cut_width(day, 1)), scale = 'width', fill = '#A40019', color = '#A40019') + 
        scale_x_continuous(breaks = -1:10) +
		theme_bw() + labs(x = 'Day', y = expression(~S[obs]))

# plot inverse simpson by day
alpha_invsimp_plot <- alpha_df %>% 
	select(group, invsimpson) %>% 
	left_join(select(meta_df, group, day), by = c('group')) %>% 
	ggplot(aes(x = day, y = invsimpson)) + 
		geom_violin(aes(group = cut_width(day, 1)), scale = 'width', fill = '#A40019', color = '#A40019') + 
        scale_x_continuous(breaks = -1:10) +
		theme_bw() + labs(x = 'Day', y = 'Inverse Simpson')

# plot theta yc by day
beta_plot <- beta_df %>% 
	inner_join(select(meta_df, group, mouse_id, day), by = c('rows' = 'group')) %>% 
	inner_join(select(meta_df, group, mouse_id, day), by = c('columns' = 'group')) %>% 
	filter(mouse_id.x == mouse_id.y, 
		day.x == -1) %>% 
		ggplot(aes(x = day.y, y = distances, group = mouse_id.x)) + 
	        scale_x_continuous(breaks = -1:10) +
	        coord_cartesian(ylim = c(0,1)) +
			geom_violin(aes(group = cut_width(day.y, 1)), scale = 'width', fill = '#A40019', color = '#A40019') + 
			theme_bw() + 
			labs(x = 'Day', y = 'Theta yc')

# save plot, top row is colonization plot, middle row are diversity plots, bottom row is temporal abundance plot
ggsave('results/figures/figure_1.jpg', plot_grid(colonization_plot, 
	plot_grid(alpha_sobs_plot, alpha_invsimp_plot, beta_plot, nrow = 1, labels = c('B', 'C', 'D')), 
	abundance_plot, ncol = 1, labels = c('A', '', 'D'), rel_heights = c(1, 1, 2)), width = 10, height = 10)
