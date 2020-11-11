##############
#
# run script to generate plots for Figure 3
#	What associates with clearance?
# 
# Nick Lesniak 04-08-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/sample.final.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#	data/mothur/sample.final.groups.ave-std.summary
#	data/mothur/sample.final.thetayc.0.03.lt.ave.dist
#	code/read.dist.R
#
##############


library(tidyverse)
library(cowplot) 
library(grid) # convert plot to grob to edit individual facet labels
library(ggpubr) # as_ggplot to convert grob to ggplot for creating figure panels
library(ggtext)  # element_markdown to allow italicized and normal text on axis ticks

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
sum_taxa_function <- 'code/sum_otu_by_taxa.R'
shared_file <- 'data/mothur/sample.final.0.03.subsample.shared'
alpha_div_file <- 'data/mothur/sample.final.groups.ave-std.summary'

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

# colors for each antibiotic
abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))
# function to background of facets
edit_facet_background <- function(plot, fills){
	#Generate the ggplot2 plot grob
	g <- grid.force(ggplotGrob(plot))
	# Get the names of grobs and their gPaths into a data.frame structure
	grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
	# Build optimal gPaths that will be later used to identify grobs and edit them
	grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
	grobs_df$gPath_full <- gsub("layout::", "",  grobs_df$gPath_full, fixed = TRUE)
	# Get the gPaths of the strip background grobs
	strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = "strip-r.*.background.*", 
		x = grobs_df$gPath_full)]
	# Get the gPaths of the strip titles
	strip_txt_gpath <- grobs_df$gPath_full[grepl(pattern = "strip-r.*titleGrob.*text.*", 
		x = grobs_df$gPath_full)]
	# Edit the grobs
	for (i in 1:length(strip_bg_gpath)){
		g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
		g <- editGrob(grob = g, gPath = strip_txt_gpath[i], gp = gpar(col = 'white'))
	}
	as_ggplot(g)
}

# get minimum relative abundance to set out limit of detection in plots
n_seqs <- median(apply(shared_df[,-1], 1, sum))
min_rel_abund <- 100 * 1/n_seqs
lod_df <- data.frame(x = 0.75, y = min_rel_abund)

# create dataframe with relative abundance subset with samples used in following analysis
meta_abund_df <- meta_df %>%
	filter(cdiff == T, time_point != 'Intermediate') %>%  # filter mice the were challenged and samples from day of challenge
	mutate(time_point = ifelse(time_point == 'Day 0', 'TOC', time_point)) %>% 
	inner_join(shared_df, by = c('group' = 'Group')) %>% 
	pivot_longer(names_to = 'OTU', values_to = 'abundance', starts_with('Otu')) %>% # convert abundances to single column
	select(group, day, abx, mouse_id, OTU, clearance, time_point, abundance) %>% 
	group_by(mouse_id, time_point) %>% 
	mutate(total = sum(abundance),
		abundance = abundance/total * 100) %>% # generate relative abundance
	group_by(abx, OTU) %>% 
	ungroup

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
		time_point = factor(time_point, levels = c('Initial', 'TOC', 'End'),
			labels = c('Initial', 'Time of challenge', 'End of experiment'))) 
# create label df to eliminate over plotting of labels
colon_clear_otu_label <- pval_diff_colon_clear_df %>% 
		select(abx, order, tax_otu_label) %>% 
		mutate(tax_otu_label = gsub(' \\(', '* \\(', tax_otu_label),
			tax_otu_label = paste0('*', tax_otu_label),
			tax_otu_label = gsub('_', ' ', tax_otu_label)) %>% 
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
	ggplot(aes(x = -order, color = time_point)) + 
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
		geom_segment(aes(y = Cleared, yend = Colonized, xend = -order)) +
		geom_point(aes(y = (abundance) + 0.04, shape = clearance), 
			position = position_dodge(width = .7), alpha = 0.2) + 
		geom_point(aes(y = Cleared), shape = 1, stroke = 1, size = 3) + 
		geom_point(aes(y = Colonized), shape = 16, size = 3) + 
		scale_shape_manual(values = c(1,16), breaks = c('Cleared', 'Colonized')) + 
		scale_color_manual(values = c('green4', 'blue3', 'red3'),
			breaks = c('Initial', 'Time of challenge', 'End of experiment')) + 
		# plot layout
		scale_y_log10(limits = c(0.04,100),
	   		breaks = c(0.01, 0.1, 1, 10, 100),
	   		labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) + 
		coord_flip() + theme_bw() + 
		labs(x = NULL, y = 'Relative Abundance (%)', shape = 'Outcome', color = 'Time Point') + 
		theme(panel.grid.minor.x = element_blank(),
			legend.position = 'bottom', 
			legend.background = element_rect(colour = 'black'),
			panel.spacing.y = unit(1.5, 'lines'),
			text = element_text(size = 10), 
			axis.text.y = element_markdown(), axis.text.x = element_markdown()) + 
		guides(shape = guide_legend(override.aes = list(alpha = 1)),
			color = guide_legend(override.aes = list(linetype = 0))) +
		facet_grid(abx~time_point, scales = 'free', space = 'free') 
# edit color of facet label background
fills <- c(pull(filter(abx_color, abx == 'Cefoperazone'), color),
		pull(filter(abx_color, abx == 'Streptomycin'), color))
diff_abund_clear_colon_plot <- edit_facet_background(diff_abund_clear_colon_plot, fills)
# add labels
diff_abund_clear_colon_plot <- plot_grid(
	plot_grid(NULL, NULL, labels = c('A', 'B'), ncol = 1, rel_heights = c(8,5)),
	diff_abund_clear_colon_plot, nrow = 1, rel_widths = c(1,19))

#### Compare the changes with in a mouse between time points , split by end status #####


# get OTUs significantly different between day of challenge and at the end of the experiment for mice that clear
pval_diff_cleared_df <- meta_abund_df %>% 
	filter(clearance %in% c('Cleared', 'Colonized')) %>% # only mice that colonization clears
	select(-day, -group, -total) %>% 
	pivot_wider(names_from = 'time_point', values_from = 'abundance') %>% 
	group_by(abx, OTU, clearance) %>% # compare OTUs across timepoints within each antibiotic
	mutate(median_TOC = median(TOC, na.rm = T), median_end = median(End, na.rm = T), median_in = median(Initial, na.rm = T)) %>%
	filter(median_TOC > 0.5 | median_end > 0.5 | median_in > 0.5)  %>% # select otus with median relative abundance > 0.5%
	nest() %>% 
	mutate(TOC_End = map(.x = data, .f = ~wilcox.test(.x$TOC, .x$End)$p.value), # compare time of challenge to end point
		Initial_TOC = map(.x = data, .f = ~wilcox.test(.x$Initial, .x$TOC)$p.value)) %>% # compare time of challenge to initial
		# did not do paired comparison because strep is missing 4 initial samples and one TOC sample resulting in a loss of all significant results
		# only comparison gained wit paired is Clindamycin OTU 15 initial VS TOC, which is initial present (~3.5%) in 2 cages, 7 mice, 
		# decreases in all cages (~0.1% RA), and recovers in one cage but not the other
	unnest(TOC_End, Initial_TOC) %>% 
	pivot_longer(names_to = 'comparison', values_to = 'pvalue', c(TOC_End, Initial_TOC)) %>% # combine all pvalues into one column
	group_by(abx) %>% 
	mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% # adjust pvalue
	filter(pvalue < 0.05) %>% # filter only those below 0.05 after correction
	unnest(data) %>% 
	pivot_longer(names_to = 'time_point', values_to = 'abundance', c(Initial, TOC, End)) %>% # unnest abundance
	filter(!is.na(abundance)) %>% 
	filter(comparison == 'TOC_End' & time_point %in% c('TOC', 'End') | 
		comparison == 'Initial_TOC' & time_point %in% c('TOC', 'Initial')) %>% # only keep abundances that match comparison 
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
			tax_otu_label = gsub('_', ' ', tax_otu_label),
			tax_otu_label = paste0('*', tax_otu_label),
			tax_otu_label = ifelse(both == T, paste0('**', tax_otu_label, '**'), 
				tax_otu_label)) %>% 
		unique

	point_shape <- ifelse(end_status == 'Cleared', 1, 16)
	point_stroke <- ifelse(end_status == 'Cleared', 1, 2)

	output_plot <- plot_df %>% 
		mutate(Outcome = as.character(point_shape)) %>% 
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
			geom_segment(aes(y = Initial, yend = TOC, xend = -order), 
					arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.5) + 
			geom_segment(aes(y = TOC, yend = End, xend = -order), 
					arrow = arrow(type = 'closed', angle = 10), color = 'black', size = 0.25) + 
			geom_point(aes(y = (abundance) + 0.04, , shape = Outcome), 
				position = position_dodge(width = .7), alpha = 0.3) + 
			geom_point(aes(y = Initial), color = 'green4', size = 3, 
				shape = point_shape, stroke = point_stroke) + 
			geom_point(aes(y = TOC), color = 'blue3', size = 3, 
				shape = point_shape, stroke = point_stroke) + 
			geom_point(aes(y = End), color = 'red3', size = 3, 
				shape = point_shape, stroke = point_stroke) + 
			scale_color_manual(values = c('green4', 'blue3', 'red3'),
				breaks = c('Initial', 'TOC', 'End')) + 
			# plot layout
			scale_y_log10(limits = c(0.04,100),
		   		breaks = c(0.01, 0.1, 1, 10, 100),
		   		labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) + 
			coord_flip() + theme_bw() +  
			labs(x = NULL, y = 'Relative Abundance (%)', color = 'Time Point') + 
			theme(panel.grid.minor = element_blank(),
				legend.position = 'none', 
				text = element_text(size = 10),
				axis.text.y = element_markdown(), axis.text.x = element_markdown())
}

###############################################################################
# plot difference between time points of mice that cleared C diff
###############################################################################
diff_abund_cleared_plot <- plot_temporal_diff_by_clearance(end_status = 'Cleared', 
	antibiotics = c('Clindamycin', 'Streptomycin'), lod_label_df = main_lod_label_df) + 
	facet_grid(abx~comparison, scales = 'free_y', space = 'free',
		labeller = labeller(comparison = c(Initial_TOC = "Initial vs Time of challenge", TOC_End = "Time of challenge vs End of experiment"))) + 
	scale_shape_manual(values = c(1,16), limits = c('Cleared', 'Colonized')) + 
	theme(panel.spacing.y = unit(1, 'lines'),
		legend.position = 'right', 
		legend.background = element_rect(colour = 'black')) + 
	guides(colour = guide_legend(override.aes = list(alpha = 1)),
		shape = guide_legend(override.aes = list(alpha = 1)))
# edit colors of facet label backgrounds
fills <- c(pull(filter(abx_color, abx == 'Clindamycin'), color),
		pull(filter(abx_color, abx == 'Streptomycin'), color))
diff_abund_cleared_plot <- edit_facet_background(diff_abund_cleared_plot, fills)
# attach labels to this part of figure
diff_abund_cleared_plot <- plot_grid(
	plot_grid(NULL, NULL, labels = c('A', 'B'), ncol = 1, rel_heights = c(9,5)),
	diff_abund_cleared_plot, nrow = 1, rel_widths = c(1,19))

###############################################################################
# plot difference between time points of mice that remained colonized by C diff
###############################################################################
diff_abund_colon_plot <- plot_temporal_diff_by_clearance(end_status = 'Colonized',
	antibiotics = c('Cefoperazone', 'Streptomycin'), lod_label_df = main_lod_label_df) + 
	facet_grid(abx~comparison, scales = 'free_y', space = 'free',
		labeller = labeller(comparison = c(Initial_TOC = "Initial vs Time of challenge", TOC_End = "Time of challenge vs End of experiment"))) + 
	theme(panel.spacing.y = unit(1, 'lines'))
# edit colors of facel label backgrounds
fills <- c(pull(filter(abx_color, abx == 'Cefoperazone'), color),
		pull(filter(abx_color, abx == 'Streptomycin'), color))
diff_abund_colon_plot <- edit_facet_background(diff_abund_colon_plot, fills)

# attach labels to this part of figure
diff_abund_colon_plot <- plot_grid(
	plot_grid(NULL, NULL, labels = c('C', 'D'), ncol = 1, rel_heights = c(8,5)),
	diff_abund_colon_plot, nrow = 1, rel_widths = c(1,19))

# plot single difference between time points of mice that cleared C diff for cef
cef_cleared_supp_plot <- plot_temporal_diff_by_clearance(end_status = 'Cleared', 
	antibiotics = c('Cefoperazone'), lod_label_df = cef_lod_label_df) + 
	facet_grid(abx~comparison, scales = 'free_y', space = 'free',
		labeller = labeller(comparison = c(Initial_TOC = "Initial vs Time of challenge", TOC_End = "Time of challenge vs End of experiment"))) +
	theme(legend.position = 'bottom',
		legend.background = element_rect(colour = 'black')) + 
	guides(shape = 'none',
		colour = guide_legend(override.aes = list(alpha = 1)))
fills <- c(pull(filter(abx_color, abx == 'Cefoperazone'), color))
cef_cleared_supp_plot <- edit_facet_background(cef_cleared_supp_plot, fills)


ggsave('results/figures/figure_3.jpg', 
	diff_abund_clear_colon_plot,
	width = 10, height = 10, units = 'in')

ggsave('results/figures/figure_4.jpg', 
	plot_grid(diff_abund_cleared_plot, diff_abund_colon_plot, nrow = 1, rel_widths = c(8,7)),
	width = 20, height = 10, units = 'in')

ggsave('results/figures/figure_S2.jpg', cef_cleared_supp_plot, 
	width = 6, height = 4, units = 'in')
