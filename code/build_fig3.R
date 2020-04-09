##############
#
# run script to generate plots for Figure 3
#	What assocaites with clearance?
# 
# Nick Lesniak 04-08-2020
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

# get minimum relative abundance to set out limit of detection in plots
n_seqs <- median(apply(shared_df[,-1], 1, sum))
min_rel_abund <- 100 * 1/n_seqs
lod_df <- data.frame(x = 0.75, y = min_rel_abund)

# create dataframe with relative abundance subset with samples used in following analysis
meta_abund_df <- meta_df %>%
	filter(cdiff == T, time_point != 'Intermediate') %>%  # filter mice the were challenged and samples from day of infection
	mutate(time_point = ifelse(time_point == 'Day 0', 'TOI', time_point)) %>% 
	inner_join(shared_df, by = c('group' = 'Group')) %>% 
	pivot_longer(names_to = 'OTU', values_to = 'abundance', starts_with('Otu')) %>% # convert abundances to single column
	select(day, abx, mouse_id, OTU, clearance, time_point, abundance) %>% 
	group_by(mouse_id, time_point) %>% 
	mutate(total = sum(abundance),
		abundance = abundance/total * 100) %>% # generate relative abundance
	group_by(abx, OTU) %>% 
	ungroup

# get OTUs significantly different between colonized and cleared by antibiotic and time point
pval_diff_colon_clear_df <- meta_abund_df %>% 
	filter(abx != 'Clindamycin', 
		clearance %in% c('Colonized', 'Cleared')) %>%  # Only compare mice that are either colonized or cleared
	pivot_wider(names_from = 'clearance', values_from = 'abundance') %>% 
	group_by(abx, time_point, OTU) %>% # compare each OTU abundance within each antibiotic in each time point 
	mutate(median_col = median(Colonized, na.rm = T), median_cle = median(Cleared, na.rm = T)) %>%
	filter(median_col > 0.5 | median_cle > 0.5)  %>% # select otus with median relative abundance > 0.5%
	nest() %>% 
	mutate(pvalue = map(.x = data, .f = ~wilcox.test(.x$Cleared, .x$Colonized)$p.value)) %>% # compare cleared vs colonized
	unnest(pvalue) %>% 
	group_by(abx, time_point, OTU) %>% 
	mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% # correct p values
	filter(pvalue < 0.05) %>% # select only those above 0.05 after pvalue correction
	unnest(data) %>% 
	pivot_longer(names_to = 'clearance', values_to = 'abundance', c(Colonized, Cleared)) %>% # unnest and recreate abundance column
	filter(!is.na(abundance)) %>% 
	left_join(select(tax_df, OTU, tax_otu_label), by = 'OTU') %>% # add otu labels
	group_by(abx, time_point, OTU) %>% 
	mutate(order = mean(abundance)) %>% # find mean abundance of comparsions to set order in plot
	ungroup
# create median df for plot medain data and line of difference 
colon_clear_median_df <- pval_diff_colon_clear_df %>% 
	group_by(abx, time_point, tax_otu_label, clearance) %>% 
	summarise(median = median(abundance) + 0.04) %>% 
	pivot_wider(names_from = clearance, values_from = median) %>% 
	ungroup

# get OTUs significantly different between day of infection and at the end of the experiment for mice that clear
pval_diff_cleared_df <- meta_abund_df %>% 
	filter(clearance == 'Cleared') %>% # only mice that colonization clears
	pivot_wider(names_from = 'time_point', values_from = 'abundance') %>% 
	group_by(abx, OTU) %>% # compare OTUs across timepoints within each antibiotic
	mutate(median_TOI = median(TOI, na.rm = T), median_end = median(End, na.rm = T), median_in = median(Initial, na.rm = T)) %>%
	filter(median_TOI > 0.5 | median_end > 0.5 | median_in > 0.5)  %>% # select otus with median relative abundance > 0.5%
	group_by(abx, OTU) %>% # compare OTUs across timepoints within each antibiotic
	nest() %>% 
	mutate(TOI_End = map(.x = data, .f = ~wilcox.test(.x$TOI, .x$End)$p.value), # compare time of infection to end point
		Initial_TOI = map(.x = data, .f = ~wilcox.test(.x$Initial, .x$TOI)$p.value)) %>% # compare time of infection to initial
	unnest(TOI_End, Initial_TOI) %>% 
	pivot_longer(names_to = 'comparison', values_to = 'pvalue', c(TOI_End, Initial_TOI)) %>% # combine all pvalues into one column
	group_by(abx, OTU, comparison) %>% 
	mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% # adjust pvalue
	filter(pvalue < 0.05) %>% # filter only those below 0.05 after correction
	unnest(data) %>% 
	pivot_longer(names_to = 'time_point', values_to = 'abundance', c(Initial, TOI, End)) %>% # unnest abundance
	filter(!is.na(abundance)) %>% 
	filter(comparison == 'TOI_End' & time_point %in% c('TOI', 'End') | 
		comparison == 'Initial_TOI' & time_point %in% c('TOI', 'Initial')) %>% # only keep abundances that match comparison 
	left_join(select(tax_df, OTU, tax_otu_label), by = 'OTU') %>% # add otu labels
	group_by(abx, comparison, OTU) %>% 
	mutate(order = mean(abundance)) %>% # find mean abundance of comparsions to set order in plot
	ungroup

cleared_median_df <- pval_diff_cleared_df %>% 
	group_by(abx, comparison, tax_otu_label, time_point) %>% 
	summarise(median = median(abundance) + 0.04) %>% 
	pivot_wider(names_from = time_point, values_from = median) %>% 
	ungroup

diff_abund_cleared_plot <- pval_diff_cleared_df %>% 
	full_join(cleared_median_df, by = c('abx', 'comparison', 'tax_otu_label')) %>% 
	mutate(tax_otu_label = gsub('_unclassified', '', tax_otu_label)) %>% 
	ggplot(aes(x = reorder(tax_otu_label, -order), color = time_point)) + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
			linetype = 'solid', color = 'black') + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
			linetype = 'dashed', color = 'white') + 
		geom_segment(aes(y = Initial, yend = TOI, 
				xend = reorder(tax_otu_label, -order)), 
				arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.5) + 
		geom_segment(aes(y = TOI, yend = End,
				xend = reorder(tax_otu_label, -order)), 
				arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.25) + 
		geom_point(aes(y = (abundance) + 0.04), 
			position = position_dodge(width = .7), alpha = 0.2) + 
		geom_point(aes(y = Initial), color = 'green4', size = 3) + 
		geom_point(aes(y = TOI), color = 'blue3', size = 3) + 
		geom_point(aes(y = End), color = 'red3', size = 3) +
		scale_y_log10(
	   		breaks = scales::trans_breaks("log10", function(x) 10^x),
	   		labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
		coord_flip() + theme_bw() +  
		labs(x = NULL, y = 'Relative Abundance (%)', color = 'Time Point',
			title = 'Differences in time of communities able to clear colonization',
			caption = 'Only significant comparisons plotted (p < 0.05 after Benjamini & Hochberg correction)') + 
		theme(legend.position = c(0.925, 0.925), 
			legend.background = element_rect(color = "black"),
			legend.title = element_text(size = 8),
			legend.text = element_text(size = 6)) + 
		geom_label(data = lod_df, aes(x = x, y = y), label = "LOD", 
			fill = "white", color = 'black', label.size = NA, inherit.aes = FALSE) + 
		facet_grid(abx~comparison, scales = 'free_y', space = 'free',
			labeller = labeller(comparison = c(Initial_TOI = "Initial vs Time of Infection", TOI_End = "Time of Infection vs End of experiment"))) + 
		theme(text = element_text(size = 10)) + 
		guides(colour = guide_legend(override.aes = list(alpha = 1)))

diff_abund_clear_colon_plot <- pval_diff_colon_clear_df %>% 
	full_join(colon_clear_median_df, by = c('abx', 'tax_otu_label', 'time_point')) %>% 
	mutate(tax_otu_label = gsub('_unclassified', '', tax_otu_label),
		time_point = factor(time_point, levels = c('Initial', 'TOI', 'End'),
			labels = c('Initial', 'Time of infection', 'End of experiment'))) %>% 
	ggplot(aes(x = reorder(tax_otu_label, -order), color = clearance)) + 
		geom_segment(aes(y = Cleared, yend = Colonized, 
			xend = reorder(tax_otu_label, -order)), color = 'black') +
		geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
			linetype = 'solid', color = 'black') + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
			linetype = 'dashed', color = 'white') + 
		geom_point(aes(y = (abundance) + 0.04), 
			position = position_dodge(width = .7), alpha = 0.2) + 
		geom_point(aes(y = Cleared), color = 'red3', size = 3) + 
		geom_point(aes(y = Colonized), color = 'cyan4', size = 3) + 
		scale_y_log10(
	   		breaks = c(0.01, 0.1, 1, 10, 100),
	   		labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
		coord_flip() + theme_bw() +  
		labs(x = NULL, y = 'Relative Abundance (%)', color = 'End Status',
			title = 'Difference between communities able and unable clear colonization',
			caption = 'Only significant comparisons plotted (p < 0.05 after Benjamini & Hochberg correction)') + 
		theme(legend.position = c(0.15, 0.08), 
			legend.background = element_rect(color = "black"),
			legend.title = element_text(size = 8),
			legend.text = element_text(size = 6)) + 
		geom_label(data = lod_df, aes(x = x, y = y), label = "LOD", 
			fill = "white", color = 'black', label.size = NA, inherit.aes = FALSE) + 
		facet_grid(abx~time_point, scales = 'free_y', space = 'free') +
		theme(text = element_text(size = 10)) + 
		guides(colour = guide_legend(override.aes = list(alpha = 1))) 

ggsave('results/figures/figure_3_diff_abund_plot.jpg', plot_grid(diff_abund_clear_colon_plot, diff_abund_cleared_plot, labels = c('A', 'B')), 
	width = 15, height = 10, units = 'in')



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
cef_abundance <- plot_day0_abundance('Cefoperazone', 12)
strep_abundance <- plot_day0_abundance('Streptomycin', 12)

plot_grid(strep_cfu_plot, cef_cfu_plot)