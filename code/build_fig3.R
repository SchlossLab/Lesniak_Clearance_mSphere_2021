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
library(ggtext)  # remotes::install_github("wilkelab/ggtext")

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

###############################################################################
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
	mutate(time_point = ifelse(time_point == 'Day 0', 'TOI', time_point)) %>% 
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
		c_time_point == 'TOI' & r_time_point == 'Initial') %>% 
	mutate(comparison = 'i_t')
TOI_intra_TOI <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx == r_abx,
		c_time_point == 'TOI' & r_time_point == 'TOI') %>% 
	mutate(comparison = 't_t')
TOI_inter_TOI <- meta_beta_df %>% 
	filter(r_mouse_id != c_mouse_id,
		c_abx != r_abx,
		c_time_point == 'TOI' & r_time_point == 'TOI') %>% 
	mutate(comparison = 'tXt')
end_toi <- meta_beta_df %>% 
	filter(r_mouse_id == c_mouse_id,
		c_time_point == 'End' & r_time_point == 'TOI') %>% 
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
			'End\nvs\nInitial')) %>% 
	mutate(c_abx = factor(c_abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>% 
	ggplot(aes(x = comparison, y = distances, color = as.factor(c_time_point))) + 
		coord_cartesian(ylim = c(0,1)) +
		geom_boxplot(aes(fill = c_clearance)) + 
		scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) +
		facet_wrap(.~c_abx) + 
		theme_bw() + 
		labs(x = NULL, y = expression(theta[YC]), fill = 'Outcome', color = 'Time Point') + 
		scale_color_manual(values = c('green4', 'blue3', 'red3'),
			breaks = c('Initial', 'TOI', 'End')) + 
		theme(legend.position = 'bottom',
			legend.key.size = unit(0.2, 'in'),
			legend.background = element_rect(color = "black"),
			panel.spacing = unit(c(3,3),'lines'))
# add labels to plots for figure
beta_plot <- plot_grid(
	plot_grid(NULL, NULL, NULL, labels = c('A', 'B', 'C'), nrow = 1), 
	beta_plot, ncol = 1, rel_heights = c(1,19))

beta_supp_plot <- beta_div_df %>% 
	filter(comparison %in% c('End\nvs\nintra\nEnd', 
			'End\nvs\ninter\nEnd')) %>% 
	mutate(c_abx = factor(c_abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin')),
		comparison = case_when(comparison == 'End\nvs\nintra\nEnd' ~ 'Within\nAntibiotic', 
			comparison == 'End\nvs\ninter\nEnd' ~ 'Across\nAntibiotic'),
		comparison = factor(comparison, levels = c('Within\nAntibiotic', 'Across\nAntibiotic'))) %>% 
	ggplot(aes(x = comparison, y = distances, color = c_abx)) + 
		coord_cartesian(ylim = c(0,1)) +
		geom_boxplot(aes(fill = c_clearance)) + 
		scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) +
		facet_grid(.~c_abx) + 
		theme_bw() + 
		labs(x = NULL, y = expression(theta[YC]), fill = 'Outcome') + 
		scale_color_manual(breaks = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
			values = c('#D37A1F', '#3A9CBC', '#A40019')) + 
		guides(color = 'none') + 
		theme(legend.position = 'bottom',
			legend.key.size = unit(0.2, 'in'),
			legend.background = element_rect(color = "black"))


###############################################################################
#   What OTUs are associated with clearance?
####

### Compare differences between communities that clear to those that cannot ###

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
	summarise(order = mean(abundance)) %>% 
	arrange(abx, order) %>% 
	ungroup %>% 
	mutate(order = row_number()) %>% # use row number inorder to maintain order within facets
	inner_join(pval_diff_colon_clear_df, by = c('abx', 'OTU'))
# create median df for plot medain data and line of difference 
colon_clear_median_df <- pval_diff_colon_clear_df %>% 
	group_by(abx, time_point, tax_otu_label, clearance) %>% 
	summarise(median = median(abundance) + 0.04) %>% 
	pivot_wider(names_from = clearance, values_from = median) %>% 
	ungroup
# create plot df with factored timepoints to correct placement
pval_diff_colon_clear_df <- pval_diff_colon_clear_df %>% 
	full_join(colon_clear_median_df, by = c('abx', 'tax_otu_label', 'time_point')) %>% 
	mutate(tax_otu_label = gsub('_unclassified', '', tax_otu_label),
		time_point = factor(time_point, levels = c('Initial', 'TOI', 'End'),
			labels = c('Initial', 'Time of infection', 'End of experiment'))) 
# create label df to eliminate over plotting of labels
colon_clear_otu_label <- pval_diff_colon_clear_df %>% 
		select(abx, order, tax_otu_label) %>% 
		mutate(tax_otu_label = gsub(' \\(', '* \\(', tax_otu_label),
			tax_otu_label = paste0('*', tax_otu_label)) %>% 
		unique

# create df to plot LOD on one set of graphs instead of all
lod_label_df <- pval_diff_colon_clear_df %>% 
	filter(abx == 'Streptomycin') %>% 
	filter(order == max(order)) %>% 
	select(abx, tax_otu_label, order) %>% 
	unique %>% 
	inner_join(distinct(select(pval_diff_colon_clear_df, abx, time_point)),
		by = c('abx')) %>% 
	mutate(y = min_rel_abund, fill = 'white', color = 'black')
# plot difference between colonized and cleared communities by abx/time point
diff_abund_clear_colon_plot <- pval_diff_colon_clear_df %>%
	ggplot(aes(x = -order)) + 
		# specify otu labels by row number
		scale_x_continuous(breaks = -colon_clear_otu_label$order, 
			labels = colon_clear_otu_label$tax_otu_label, expand = c(0,0)) + 
		# limit of detection
		geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
			linetype = 'solid', color = 'black') + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
			linetype = 'dashed', color = 'white') + 
		geom_label(data = lod_label_df, aes(y = y), label = 'LOD', color = 'white') + 
		geom_text(data = lod_label_df, aes(y = y), label = 'LOD', color = 'black') + 
		#  barbell geom
		geom_segment(aes(y = Cleared, yend = Colonized, xend = -order), color = 'black') +
		geom_point(aes(y = (abundance) + 0.04, shape = clearance), 
			position = position_dodge(width = .7), alpha = 0.2) + 
		geom_point(aes(y = Cleared), shape = 1, stroke = 1, size = 3) + 
		geom_point(aes(y = Colonized), shape = 16, size = 3) + 
		scale_shape_manual(values = c(1,16), breaks = c('Cleared', 'Colonized')) + 
		# plot layout
		scale_y_log10(limits = c(0.04,100),
	   		breaks = c(0.01, 0.1, 1, 10, 100),
	   		labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) + 
		coord_flip() + theme_bw() + 
		labs(x = NULL, y = 'Relative Abundance (%)', shape = 'Outcome') + 
		theme(panel.grid.minor.x = element_blank(),
			legend.position = 'none', 
			panel.spacing.y = unit(3, 'lines'),
			text = element_text(size = 10), 
			axis.text.y = element_markdown(), axis.text.x = element_markdown()) + 
		facet_grid(abx~time_point, scales = 'free', space = 'free') 
diff_abund_clear_colon_plot <- plot_grid(
	plot_grid(NULL, NULL, labels = c('D', 'E'), ncol = 1, rel_heights = c(9,3.5)),
	diff_abund_clear_colon_plot, nrow = 1, rel_widths = c(1,19))

#### Compare the changes with in a mouse between time points , split by end status #####


# get OTUs significantly different between day of infection and at the end of the experiment for mice that clear
pval_diff_cleared_df <- meta_abund_df %>% 
	filter(clearance %in% c('Cleared', 'Colonized')) %>% # only mice that colonization clears
	select(-day, -group, -total) %>% 
	pivot_wider(names_from = 'time_point', values_from = 'abundance') %>% 
	group_by(abx, OTU, clearance) %>% # compare OTUs across timepoints within each antibiotic
	mutate(median_TOI = median(TOI, na.rm = T), median_end = median(End, na.rm = T), median_in = median(Initial, na.rm = T)) %>%
	filter(median_TOI > 0.5 | median_end > 0.5 | median_in > 0.5)  %>% # select otus with median relative abundance > 0.5%
	nest() %>% 
	mutate(TOI_End = map(.x = data, .f = ~wilcox.test(.x$TOI, .x$End)$p.value), # compare time of infection to end point
		Initial_TOI = map(.x = data, .f = ~wilcox.test(.x$Initial, .x$TOI)$p.value)) %>% # compare time of infection to initial
		# did not do paired comparison because strep is missing 4 initial samples and one TOI sample resulting in a loss of all significant results
		# only comparison gained wit paired is Clindamycin OTU 15 initial VS TOI, which is initial present (~3.5%) in 2 cages, 7 mice, 
		# decreases in all cages (~0.1% RA), and recovers in one cage but not the other
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
	left_join(select(tax_df, OTU, tax_otu_label), by = 'OTU') # add otu labels
# find mean abundance of comparsions to set order in plot
pval_diff_cleared_df <- pval_diff_cleared_df %>% 
	group_by(abx, clearance, OTU) %>% 
	summarise(order = mean(abundance)) %>% 
	arrange(clearance, abx, order) %>% 
	ungroup %>% 
	mutate(order = row_number()) %>% # use absolute row number inorder to maintain order within facets
	inner_join(pval_diff_cleared_df, by = c('abx', 'clearance', 'OTU'))
# create median df for summary barbell geom
cleared_median_df <- pval_diff_cleared_df %>% 
	group_by(abx, clearance, comparison, tax_otu_label, time_point) %>% 
	summarise(median = median(abundance) + 0.04) %>% 
	pivot_wider(names_from = time_point, values_from = median) %>% 
	ungroup
# create df to plot LOD on one set of graphs instead of all
main_lod_label_df <- pval_diff_cleared_df %>% 
	filter(abx == 'Streptomycin') %>% 
	group_by(clearance) %>% 
	filter(order == max(order)) %>% ungroup %>% 
	select(abx, tax_otu_label, clearance, order) %>% 
	unique %>% 
	inner_join(distinct(select(pval_diff_cleared_df, abx, comparison)),
		by = c('abx')) %>% 
	mutate(y = min_rel_abund, fill = 'white', color = 'black', time_point = NA)

cef_lod_label_df <- pval_diff_cleared_df %>% 
	filter(abx == 'Cefoperazone') %>% 
	group_by(clearance) %>% 
	filter(order == max(order)) %>% ungroup %>% 
	select(abx, tax_otu_label, clearance, order) %>% 
	unique %>% 
	mutate(order = order + 0.25) %>% 
	inner_join(distinct(select(pval_diff_cleared_df, abx, comparison)),
		by = c('abx')) %>% 
	mutate(y = min_rel_abund, fill = 'white', color = 'black', time_point = NA)

# create column noting which OTUs are in both comparisons
pval_diff_cleared_df <- pval_diff_cleared_df %>% 
	select(abx, clearance, OTU) %>% 
	distinct %>% 
	group_by(abx) %>% 
	count(OTU) %>% 
	mutate(both = n > 1) %>% 
	full_join(pval_diff_cleared_df, by = c('abx', 'OTU')) %>% 
	ungroup

# plot differences between time points by clearance/time/abx
plot_temporal_diff_by_clearance <- function(end_status, antibiotics, lod_label_df){
	plot_df <- pval_diff_cleared_df %>% 
		full_join(cleared_median_df, by = c('abx', 'clearance', 'comparison', 'tax_otu_label')) %>% 
		filter(clearance == end_status,
			abx %in% antibiotics) %>% 
		mutate(tax_otu_label = gsub('_unclassified', '', tax_otu_label))
	# create label df to eliminate over plotting of labels and bold OTUs in both comparisons
	otu_label <- plot_df %>% 
		select(abx, comparison, order, tax_otu_label, both) %>% 
		mutate(tax_otu_label = gsub(' \\(', '* \\(', tax_otu_label),
			tax_otu_label = paste0('*', tax_otu_label),
			tax_otu_label = ifelse(both == T, paste0('**', tax_otu_label, '**'), 
				tax_otu_label)) %>% 
		unique

	point_shape <- ifelse(end_status == 'Cleared', 1, 16)
	point_stroke <- ifelse(end_status == 'Cleared', 1, 2)

	output_plot <- plot_df %>% 
		ggplot(aes(x = -order, color = time_point)) + 
			# specify labels for row numbers
			scale_x_continuous(breaks = -otu_label$order, 
				labels = otu_label$tax_otu_label, expand = c(0,0)) + 
			# limit of detection line
			geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
				linetype = 'solid', color = 'black') + 
			geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
				linetype = 'dashed', color = 'white') + 
			geom_label(data = filter(lod_label_df, clearance == end_status), 
				aes(y = y), label = 'LOD', color = 'white') + 
			geom_text(data = filter(lod_label_df, clearance == end_status), 
				aes(y = y), label = 'LOD', color = 'black') + 
			# points with arrows indicating direction of change
			geom_segment(aes(y = Initial, yend = TOI, xend = -order), 
					arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.5) + 
			geom_segment(aes(y = TOI, yend = End, xend = -order), 
					arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.25) + 
			geom_point(aes(y = (abundance) + 0.04), 
				position = position_dodge(width = .7), alpha = 0.3) + 
			geom_point(aes(y = Initial), color = 'green4', size = 3, 
				shape = point_shape, stroke = point_stroke) + 
			geom_point(aes(y = TOI), color = 'blue3', size = 3, 
				shape = point_shape, stroke = point_stroke) + 
			geom_point(aes(y = End), color = 'red3', size = 3, 
				shape = point_shape, stroke = point_stroke) + 
			scale_color_manual(values = c('green4', 'blue3', 'red3'),
				breaks = c('Initial', 'TOI', 'End')) + 
			# plot layout
			scale_y_log10(limits = c(0.04,100),
		   		breaks = c(0.01, 0.1, 1, 10, 100),
		   		labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) + 
			coord_flip() + theme_bw() +  
			labs(x = NULL, y = 'Relative Abundance (%)') + 
			theme(panel.grid.minor = element_blank(),
				legend.position = 'none', 
				text = element_text(size = 10),
				axis.text.y = element_markdown(), axis.text.x = element_markdown())
}
# plot difference between time points of mice that cleared C diff
diff_abund_cleared_plot <- plot_temporal_diff_by_clearance(end_status = 'Cleared', 
	antibiotics = c('Clindamycin', 'Streptomycin'), lod_label_df = main_lod_label_df) + 
	facet_grid(abx~comparison, scales = 'free_y', space = 'free',
		labeller = labeller(comparison = c(Initial_TOI = "Initial vs Time of Infection", TOI_End = "Time of Infection vs End of experiment"))) + 
	theme(panel.spacing.y = unit(3, 'lines'))
# attach labels to this part of figure
diff_abund_cleared_plot <- plot_grid(
	plot_grid(NULL, NULL, labels = c('F', 'G'), ncol = 1, rel_heights = c(8,4)),
	diff_abund_cleared_plot, nrow = 1, rel_widths = c(1,19))
# plot single difference between time points of mice that cleared C diff for cef
cef_cleared_supp_plot <- plot_temporal_diff_by_clearance(end_status = 'Cleared', 
	antibiotics = c('Cefoperazone'), lod_label_df = cef_lod_label_df) + 
	facet_grid(.~comparison, scales = 'free_y', space = 'free',
		labeller = labeller(comparison = c(Initial_TOI = "Initial vs Time of Infection", TOI_End = "Time of Infection vs End of experiment"))) 
# plot difference between time points of mice that remained colonized by C diff
diff_abund_colon_plot <- plot_temporal_diff_by_clearance(end_status = 'Colonized',
	antibiotics = c('Cefoperazone', 'Streptomycin'), lod_label_df = main_lod_label_df) + 
	facet_grid(abx~comparison, scales = 'free_y', space = 'free',
		labeller = labeller(comparison = c(Initial_TOI = "Initial vs Time of Infection", TOI_End = "Time of Infection vs End of experiment"))) + 
	theme(panel.spacing.y = unit(3, 'lines'))
# attach labels to this part of figure
diff_abund_colon_plot <- plot_grid(
	plot_grid(NULL, NULL, labels = c('H', 'I'), ncol = 1, rel_heights = c(9,4)),
	diff_abund_colon_plot, nrow = 1, rel_widths = c(1,19))


ggsave('results/figures/figure_3.jpg', 
	plot_grid(
		plot_grid(beta_plot, diff_abund_clear_colon_plot, nrow = 1),
		plot_grid(diff_abund_cleared_plot, diff_abund_colon_plot, nrow = 1),
		ncol = 1), 
	width = 25, height = 20, units = 'in')

ggsave('results/figures/figure_S1.jpg', beta_supp_plot, 
	width = 10, height = 10, units = 'in')

ggsave('results/figures/figure_S2.jpg', cef_cleared_supp_plot, 
	width = 6, height = 2, units = 'in')
