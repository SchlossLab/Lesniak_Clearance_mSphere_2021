##############
#
# run script to generate plots for Figure 3
#	What associates with clearance?
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
	select(group, day, abx, mouse_id, OTU, clearance, time_point, abundance) %>% 
	group_by(mouse_id, time_point) %>% 
	mutate(total = sum(abundance),
		abundance = abundance/total * 100) %>% # generate relative abundance
	group_by(abx, OTU) %>% 
	ungroup

####
#   How different are communities at TOI and End?
####

# plot Theta yc differences for 
#	initial as control range, 
#	TOI vs own initial, TOI vs other Initial, 
#	End vs End, End vs Initial, End vs TOI, End vs other End
beta_meta_df <- meta_df %>% 
	filter(clearance %in% c('Cleared', 'Colonized'),
		time_point != 'Intermediate',
		cdiff == TRUE) %>% 
	select(group, abx, dose, clearance, time_point, mouse_id) 

meta_beta_df <- beta_df %>% 
	inner_join(rename_all(beta_meta_df, list(~paste0('r_', .))), by = c('rows' = 'r_group')) %>% 
	inner_join(rename_all(beta_meta_df, list(~paste0('c_', .))), by = c('columns' = 'c_group')) %>% 
	filter()

initial_distances <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx == r_abx,
		r_time_point == 'Initial' & c_time_point == 'Initial') %>% 
	mutate(comparison = 'i_i',)
initial_inter_intial <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx != r_abx,
		c_time_point == 'Initial' & r_time_point == 'Initial') %>% 
	mutate(comparison = 'iXi')
TOI_intra_initial <- meta_beta_df %>% 
	filter(r_mouse_id == c_mouse_id,
		c_time_point == 'Day 0' & r_time_point == 'Initial') %>% 
	mutate(comparison = 'i_t')
TOI_intra_TOI <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx == r_abx,
		c_time_point == 'Day 0' & r_time_point == 'Day 0') %>% 
	mutate(comparison = 't_t')
TOI_inter_TOI <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx != r_abx,
		c_time_point == 'Day 0' & r_time_point == 'Day 0') %>% 
	mutate(comparison = 'tXt')
end_toi <- meta_beta_df %>% 
	filter(r_mouse_id == c_mouse_id,
		c_time_point == 'End' & r_time_point == 'Day 0') %>% 
	mutate(comparison = 'e_t')
end_intial <- meta_beta_df %>% 
	filter(r_mouse_id == c_mouse_id,
		c_time_point == 'End' & r_time_point == 'Initial') %>% 
	mutate(comparison = 'e_i')
end_intra_end <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx == r_abx,
		c_time_point == 'End' & r_time_point == 'End') %>% 
	mutate(comparison = 'e_e')
end_inter_end <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx != r_abx,
		c_time_point == 'End' & r_time_point == 'End') %>% 
	mutate(comparison = 'eXe')

beta_div_df <- bind_rows(list(initial_distances, initial_inter_intial, TOI_intra_initial, TOI_intra_TOI, 
	TOI_inter_TOI, end_toi, end_intial, end_intra_end, end_inter_end)) %>% 
	mutate(comparison = factor(comparison, 
		levels = c('i_i', 'i_t', 'e_i', 'e_t', 'iXi', 't_t', 'tXt', 'e_e', 'eXe'),
		labels = c('Initial\nvs\nInitial', 
			'Initial\nvs\nTOI', 
			'End\nvs\nInitial', 
			'End\nvs\nTOI', 
			'Initial\nvs\ninter\nInitial',
			'TOI\nvs\nintra\nTOI', 
			'TOI\nvs\ninter\nTOI',
			'End\nvs\nintra\nEnd', 
			'End\nvs\ninter\nEnd')))

beta_plot <- beta_div_df %>% 
	filter(comparison %in% c('Initial\nvs\nInitial', 
			'Initial\nvs\nTOI', 
			'End\nvs\nInitial', 
			'TOI\nvs\nintra\nTOI',  
			'End\nvs\nintra\nEnd', 
			'TOI\nvs\ninter\nTOI',
			'End\nvs\ninter\nEnd')) %>% 
	mutate(c_abx = factor(c_abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>% 
	ggplot(aes(x = comparison, y = distances, fill = as.factor(c_time_point))) + 
		coord_cartesian(ylim = c(0,1)) +
		geom_boxplot() + 
		facet_grid(c_clearance~c_abx) + 
		theme_bw() + 
		labs(x = NULL, y = 'Theta yc', fill = 'Time Point') + 
		theme(legend.position = c(0.08, 0.25),
			legend.key.size = unit(0.2, 'in'))

####
#   What OTUs are associated with clearance?
####

# get OTUs significantly different between colonized and cleared by antibiotic and time point
pval_diff_colon_clear_df <- meta_abund_df %>% 
	select(-group, -day, -total) %>% 
	filter(abx != 'Clindamycin', 
		clearance %in% c('Colonized', 'Cleared')) %>%  # Only compare mice that are either colonized or cleared
	group_by(abx, clearance, time_point, OTU) %>% # compare each OTU abundance within each antibiotic in each time point 
	mutate(median_abundance = median(abundance)) %>%
	group_by(abx, time_point, OTU) %>% 
	filter(max(median_abundance) > 0.5)  %>% # select otus with median relative abundance > 0.5% in either group
	nest() %>% 
	mutate(pvalue = map(.x = data, .f = ~ wilcox.test(
				pull( filter(.x, clearance == 'Cleared'), abundance), 
				pull( filter(.x, clearance == 'Colonized'), abundance))$p.value)) %>% # compare cleared vs colonized
	unnest(pvalue) %>% 
	group_by(abx) %>% 
	mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% # correct p values
	filter(pvalue < 0.05) %>% # select only those above 0.05 after pvalue correction
	unnest(data) %>% 
	#pivot_longer(names_to = 'clearance', values_to = 'abundance', c(Colonized, Cleared)) %>% # unnest and recreate abundance column
	#filter(!is.na(abundance)) %>% 
	left_join(select(tax_df, OTU, tax_otu_label), by = 'OTU') # add otu labels
# find mean abundance of comparsions to set order in plot
pval_diff_colon_clear_df <- pval_diff_colon_clear_df %>% 
	group_by(abx, OTU) %>% 
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
	group_by(abx) %>% 
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

ggsave('results/figures/figure_3.jpg', 
	plot_grid(plot_grid(NULL, plot_grid(beta_plot, labels = c('A')), NULL, rel_widths = c(1, 3, 1), nrow = 1), 
		plot_grid(diff_abund_clear_colon_plot, diff_abund_cleared_plot, labels = c('B', 'C'), nrow = 1),
			ncol = 1, rel_heights = c(1,2)), 
	width = 15, height = 10, units = 'in')
# cef_cleared_df <- meta_abund_df %>% 
#	filter(abx == 'Cefoperazone', clearance == 'Cleared') %>% # only mice that colonization clears
#	select(-day, -group, -clearance, -total) %>% 
#	pivot_wider(names_from = 'time_point', values_from = 'abundance') %>% 
#	group_by(OTU) %>% # compare OTUs across timepoints within each antibiotic
#	mutate(median_TOI = median(TOI, na.rm = T), median_end = median(End, na.rm = T), median_in = median(Initial, na.rm = T)) %>%
#	filter(median_TOI > 0.5 | median_end > 0.5 | median_in > 0.5)  %>% # select otus with median relative abundance > 0.5%
#	nest() %>% 
#	mutate(TOI_End = map(.x = data, .f = ~wilcox.test(.x$TOI, .x$End)$p.value), # compare time of infection to end point
#		Initial_TOI = map(.x = data, .f = ~wilcox.test(.x$Initial, .x$TOI)$p.value))  %>% 
#	unnest(TOI_End, Initial_TOI) %>% 
#	pivot_longer(names_to = 'comparison', values_to = 'pvalue', c(TOI_End, Initial_TOI)) %>% # combine all pvalues into one column
#	ungroup %>% 
#	filter(pvalue < 0.1) %>% 	
#	unnest(data) %>% 
#	pivot_longer(names_to = 'time_point', values_to = 'abundance', c(Initial, TOI, End)) %>% # unnest abundance
#	filter(!is.na(abundance)) %>% 
#	filter(comparison == 'TOI_End' & time_point %in% c('TOI', 'End') | 
#		comparison == 'Initial_TOI' & time_point %in% c('TOI', 'Initial')) %>% # only keep abundances that match comparison 
#	left_join(select(tax_df, OTU, tax_otu_label), by = 'OTU') %>% # add otu labels
#	group_by(abx, comparison, OTU) %>% 
#	mutate(order = mean(abundance)) %>% # find mean abundance of comparsions to set order in plot
#	ungroup
#cef_cleared_median_df <- cef_cleared_df %>% 
#	group_by(abx, comparison, tax_otu_label, time_point) %>% 
#	summarise(median = median(abundance) + 0.04) %>% 
#	pivot_wider(names_from = time_point, values_from = median) %>% 
#	ungroup
#cef_diff_plot <- cef_cleared_df %>%
#	full_join(cef_cleared_median_df, by = c('abx', 'comparison', 'tax_otu_label')) %>% 
#	mutate(tax_otu_label = gsub('_unclassified', '', tax_otu_label),
#		tax_otu_label = paste0(tax_otu_label, '\n(Uncorrected p-value = ', round(pvalue, 3), ')')) %>% 
#	ggplot(aes(x = reorder(tax_otu_label, -order), color = time_point)) + 
#		geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
#			linetype = 'solid', color = 'black') + 
#		geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
#			linetype = 'dashed', color = 'white') + 
#		geom_segment(aes(y = Initial, yend = TOI, 
#				xend = reorder(tax_otu_label, -order)), 
#				arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.5) + 
#		geom_segment(aes(y = TOI, yend = End,
#				xend = reorder(tax_otu_label, -order)), 
#				arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.25) + 
#		geom_point(aes(y = (abundance) + 0.04), 
#			position = position_dodge(width = .7), alpha = 0.2) + 
#		geom_point(aes(y = Initial), color = 'green4', size = 3) + 
#		geom_point(aes(y = TOI), color = 'blue3', size = 3) + 
#		geom_point(aes(y = End), color = 'red3', size = 3) +
#		scale_y_log10(
#	   		breaks = scales::trans_breaks("log10", function(x) 10^x),
#	   		labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
#		coord_flip() + theme_bw() +  
#		labs(x = NULL, y = 'Relative Abundance (%)', color = 'Time Point',
#			title = 'Temporal Differences in Cefoperazone communities able to clear colonization',
#			caption = 'No OTUs are significant after multiple comparisons correction') + 
#		theme(legend.position = c(0.925, 0.925), 
#			legend.background = element_rect(color = "black"),
#			legend.title = element_text(size = 8),
#			legend.text = element_text(size = 6)) + 
#		geom_label(data = lod_df, aes(x = x, y = y), label = "LOD", 
#			fill = "white", color = 'black', label.size = NA, inherit.aes = FALSE) + 
#		facet_grid(abx~comparison, scales = 'free_y', space = 'free',
#			labeller = labeller(comparison = c(Initial_TOI = "Initial vs Time of Infection", TOI_End = "Time of Infection vs End of experiment"))) + 
#		theme(text = element_text(size = 10)) + 
#		guides(colour = guide_legend(override.aes = list(alpha = 1)))
#ggsave('results/figures/figure_3_cef_cleared.jpg', ceff_diff_plot, width = 8, height = 6, units = 'in')
#