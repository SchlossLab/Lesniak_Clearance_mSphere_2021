##############
#
# run script to generate plots for Figure 2
#	Can other conditions naturally clear C. difficile colonization?
#
# Nick Lesniak 04-06-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/sample.final.groups.ave-std.summary
#	data/mothur/sample.final.thetayc.0.03.lt.ave.dist
#	code/read_dist.R
#
##############



library(tidyverse)
library(cowplot)
library(ggtext)  # remotes::install_github("wilkelab/ggtext")
library(ggpubr) # as_ggplot to convert grob to ggplot for creating figure panels
library(grid) # convert plot to grob to edit individual facet labels

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
alpha_div_file <- 'data/mothur/sample.final.groups.ave-std.summary'
beta_div_file <- 'data/mothur/sample.final.thetayc.0.03.lt.ave.dist'
dist_function <- 'code/read_dist.R'

abx_color <- tibble(abx = c('Clindamycin', 'Cefoperazone', 'Streptomycin'),
	color = c('#A40019', '#3A9CBC', '#D37A1F'))
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
	strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = "strip-*.background.*",
		x = grobs_df$gPath_full)]
	# Get the gPaths of the strip titles
	strip_txt_gpath <- grobs_df$gPath_full[grepl(pattern = "strip.*titleGrob.*text.*",
		x = grobs_df$gPath_full)]
	# Edit the grobs
	for (i in 1:length(strip_bg_gpath)){
		g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
		g <- editGrob(grob = g, gPath = strip_txt_gpath[i], gp = gpar(col = 'white'))
	}
	as_ggplot(g)
}

# read in data
meta_df   <- read_tsv(meta_file) %>%
	filter(clearance %in% c('Cleared', 'Colonized'),
		cdiff == T, time_point != 'Intermediate') %>%
	mutate(time_point = ifelse(time_point == 'Day 0', 'TOC', time_point),
		time_point = factor(time_point, levels = c('Initial', 'TOC', 'End'),
			labels = c('Initial', 'TOC', 'End')),
		abx = factor(abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'),
			labels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>%
	select(group, dose, day, abx, mouse_id, clearance, time_point)

alpha_df <- read_tsv(alpha_div_file) %>%
	filter(group %in% meta_df$group,
		method == 'ave') %>%
	select(group, sobs, invsimpson) %>%
	inner_join(select(meta_df, abx, group, time_point, clearance, dose), by = c('group'))


source(dist_function) # function to read in distance file and convert from triangle to dataframe
beta_df <- read_dist(beta_div_file) %>%
	filter(rows %in% meta_df$group,
		columns %in% meta_df$group)

alpha_wilcox_df <- alpha_df %>%
	mutate(subset = paste(time_point, clearance, sep = '.')) %>%
	group_by(abx) %>%
	nest() %>%
	mutate(data = map(data, function(data) nest(group_by(data, subset)))) %>%
	mutate(data = map(data, function(nested_df){
		sobs_df <- sapply(nested_df$data, function(x) sapply(nested_df$data, function(y) wilcox.test(x$sobs,y$sobs)$p.value)) %>%
			data.frame
		sobs_df[upper.tri(sobs_df)] <- NA # set all of upper triangle to 1 to eliminate duplicate comparisons
		colnames(sobs_df) <- nested_df$subset
		sobs_df <- sobs_df %>%
			mutate(row_names = nested_df$subset) %>%
			pivot_longer(col = -row_names, names_to = 'col_names', values_to = 'pvalue') %>%
			filter(pvalue != 1, !is.na(pvalue)) %>% # eliminate all self comparisons and upper triangle
			mutate(pvalue = p.adjust(pvalue, method = 'BH'),
				metric = 'sobs') # correct p values

		invsimpson_df <- sapply(nested_df$data, function(x) sapply(nested_df$data, function(y) wilcox.test(x$sobs,y$sobs)$p.value)) %>%
			data.frame
		invsimpson_df[upper.tri(invsimpson_df)] <- NA # set all of upper triangle to 1 to eliminate duplicate comparisons
		colnames(invsimpson_df) <- nested_df$subset
		invsimpson_df <- invsimpson_df %>%
			mutate(row_names = nested_df$subset) %>%
			pivot_longer(col = -row_names, names_to = 'col_names', values_to = 'pvalue') %>%
			filter(pvalue != 1, !is.na(pvalue)) %>% # eliminate all self comparisons and upper triangle
			mutate(pvalue = p.adjust(pvalue, method = 'BH'),
				metric = 'invsimpson') # correct p values

		return(bind_rows(sobs_df, invsimpson_df))
		})) %>%
	unnest(data) %>%
	ungroup

signif_label_position <- data.frame(
	abx = c(rep('Clindamycin', 3), rep(rep(c('Cefoperazone', 'Streptomycin'), each = 3), 2)),
	subset = c(rep(paste(c('Initial', 'TOC', 'End'), 'Cleared', sep = '.'), 3),
			rep(paste(c('Initial', 'TOC', 'End'), 'Colonized', sep = '.'), 2)),
	x = c(1,2,3, rep(c(0.85, 1.85, 2.85), 2), rep(c(1.15, 2.15, 3.15), 2)))

signif_label_df <- alpha_wilcox_df %>%
	filter(pvalue < 0.05) %>%
	left_join(signif_label_position, by = c('abx', 'row_names' = 'subset')) %>%
	left_join(signif_label_position, by = c('abx', 'col_names' = 'subset'),
		suffix = c('1','2')) %>%
	mutate(label = '*',
		xnote = (x1 + x2)/2,
		y1 = c(c(190, 215, 170, 190, 170, 196, 221), # cef sobs
			c(45, 51.5, 40, 45, 40, 46.5, 53), # cef inv simpson
			c(136, 115, 100, 166, 172, 100, 178, 130, 166, 121, 160), # strep sobs
			c(36.5, 30, 25, 43, 44.5, 25, 46, 35, 43, 31.5, 41.5), # strep inv simp
			c(100, 110), # clinda sobs
			c(20, 22)), # clinda inv simpson
		ynote = y1 + c(rep(.5, 7), rep(.1, 7), rep(.5, 11), rep(.1, 11), c(.5,.5,.1,.1)),
		clearance = 'NA', time_point = 'NA',
		alpha = case_when(x1 == 1.15 & x2 == 1.85 ~ 0.4,
			x2 == 1.15 & x1 == 1.85 ~ 0.4,
			x1 == 2.15 & x2 == 2.85 ~ 0.4,
			x2 == 2.15 & x1 == 2.85 ~ 0.4,
			x2 == 0.85 & x1 == 3.15 ~ 0.4,
			x1 == 0.85 & x2 == 3.15 ~ 0.4,
			x2 == 0.85 & x1 == 2.15 ~ 0.4,
			x1 == 0.85 & x2 == 2.15 ~ 0.4,
			T ~ 1))  %>%
	filter(alpha == 1)

# plot Sobs by day
alpha_sobs_plot <- alpha_df %>%
	ggplot(aes(x = time_point, y = sobs)) +
		geom_point(aes(shape = clearance, color = abx),
			position = position_jitterdodge()) +
		scale_shape_manual(values = c(1, 16), limits = c('Cleared', 'Colonized')) +
		scale_color_manual(values = abx_color$color, limits = abx_color$abx) +
		theme_bw() + labs(x = NULL, y = expression(~S[obs])) +
		theme(panel.grid.minor = element_blank(),
		  axis.title = element_text(size = 14),
		  strip.text = element_text(size = 14),
			legend.position = 'none',
			panel.spacing = unit(c(3,3),'lines')) +
		facet_wrap(.~abx) +
		geom_text(data = filter(signif_label_df, metric == 'sobs'),
			aes(x = xnote, y = ynote, label = label)) +
		geom_segment(data = filter(signif_label_df, metric == 'sobs'),
			aes(x = x1, xend = x2, y = y1, yend = y1), size = 0.25)


alpha_sobs_plot <- edit_facet_background(alpha_sobs_plot, abx_color$color)
# add labels to plots for figure
alpha_sobs_plot <- plot_grid(
	plot_grid(NULL, NULL, NULL, labels = c('A', 'B', 'C'), nrow = 1),
	alpha_sobs_plot, ncol = 1, rel_heights = c(1.5,19))


# plot inverse simpson by day
alpha_invsimpson_plot <- alpha_df %>%
	ggplot(aes(x = time_point, y = invsimpson)) +
		geom_point(aes(shape = clearance, color = abx),
			position = position_jitterdodge()) +
		scale_shape_manual(values = c(1, 16), limits = c('Cleared', 'Colonized')) +
		scale_color_manual(values = abx_color$color, limits = abx_color$abx) +
		theme_bw() + labs(x = NULL, y = 'Inverse Simpson') +
		theme(panel.grid.minor = element_blank(),
			legend.position = 'none',
			panel.spacing = unit(c(3,3),'lines'),
			axis.title = element_text(size = 14),
			strip.background = element_blank(),
			strip.text = element_blank()) +
		facet_wrap(.~abx) +
		geom_text(data = filter(signif_label_df, metric == 'invsimpson'),
			aes(x = xnote, y = ynote, label = label)) +
		geom_segment(data = filter(signif_label_df, metric == 'invsimpson'),
			aes(x = x1, xend = x2, y = y1, yend = y1), size = 0.25)

############
# Does dose contribute to high variation alpha in colonized cef mice?
############
#test_dose_toc_df <- alpha_df %>% 
#  filter(abx == 'Cefoperazone', clearance == 'Colonized', time_point == 'TOC') %>% 
#  select(sobs, invsimpson, dose) %>% 
#  pivot_longer(c(sobs, invsimpson), names_to = 'measure') %>% 
#  mutate(sample = paste(measure, dose, sep = '_')) %>% 
#  select(sample, value) %>% 
#  pivot_wider(names_from = 'sample', values_from = 'value')
#test_dose_end_df <- alpha_df %>% 
#  filter(abx == 'Cefoperazone', clearance == 'Colonized', time_point == 'End') %>% 
#  select(sobs, invsimpson, dose) %>% 
#  pivot_longer(c(sobs, invsimpson), names_to = 'measure') %>% 
#  mutate(sample = paste(measure, dose, sep = '_')) %>% 
#  select(sample, value) %>% 
#  pivot_wider(names_from = 'sample', values_from = 'value')
#pval_cef_dose <- c(wilcox.test(test_dose_toc_df$sobs_0.3[[1]], test_dose_toc_df$sobs_0.5[[1]])$p.value,
#  wilcox.test(test_dose_toc_df$invsimpson_0.3[[1]], test_dose_toc_df$invsimpson_0.5[[1]])$p.value,
#  wilcox.test(test_dose_end_df$sobs_0.3[[1]], test_dose_end_df$sobs_0.5[[1]])$p.value,
#  wilcox.test(test_dose_end_df$invsimpson_0.3[[1]], test_dose_end_df$invsimpson_0.5[[1]])$p.value)
## 0.02857143 0.20000000 0.11111111 0.41269841
#pval_cef_dose <- p.adjust(pval_cef_dose, method = 'BH')
## 0.1142857 0.2666667 0.2222222 0.4126984
#
alpha_sobs_plot_by_dose <- alpha_df %>%
  filter(abx =='Cefoperazone') %>% 
  mutate(dose = paste(dose, 'mg/mL')) %>% 
  ggplot(aes(x = time_point, y = sobs)) +
  geom_point(aes(shape = clearance),
             color = abx_color$color[2],
             position = position_jitterdodge()) +
  scale_shape_manual(values = c(1, 16), 
                     limits = c('Cleared', 'Colonized')) +
  theme_bw() + labs(x = NULL, y = expression(~S[obs])) +
  theme(panel.grid.minor = element_blank(),
        strip.background =element_rect(fill=abx_color$color[2]), 
        axis.title = element_text(size = 14),
        strip.text = element_text(colour = 'white', size = 12),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(.~dose)


alpha_invsimpson_plot_by_dose <- alpha_df %>%
  filter(abx == 'Cefoperazone') %>% 
  ggplot(aes(x = time_point, y = invsimpson)) +
  geom_point(aes(shape = clearance),
             color = abx_color$color[2],
             position = position_jitterdodge()) +
  scale_shape_manual(values = c(1, 16), 
                     limits = c('Cleared', 'Colonized')) +
  theme_bw() + labs(x = 'Day', y = 'Inverse Simpson', shape = NULL) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'bottom',
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(.~dose)



###############################################################################
#   How different are communities at TOC and End?
####

# plot Theta yc differences for
#	initial as control range,
#	TOC vs own initial, TOC vs other Initial,
#	End vs End, End vs Initial, End vs TOC, End vs other End
beta_meta_df <- meta_df %>%
	select(group, abx, dose, clearance, time_point, mouse_id)

meta_beta_df <- beta_df %>%
	inner_join(rename_all(beta_meta_df, list(~paste0('r_', .))), by = c('rows' = 'r_group')) %>%
	inner_join(rename_all(beta_meta_df, list(~paste0('c_', .))), by = c('columns' = 'c_group')) %>%
	filter()

initial_distances <- meta_beta_df %>%
	filter(r_mouse_id != c_mouse_id,
		c_abx == r_abx,
		r_time_point == 'Initial' & c_time_point == 'Initial') %>%
	mutate(comparison = 'i_i')
initial_inter_intial <- meta_beta_df %>%
	filter(r_mouse_id != c_mouse_id,
		c_abx != r_abx,
		c_time_point == 'Initial' & r_time_point == 'Initial') %>%
	mutate(comparison = 'iXi')
TOC_intra_initial <- meta_beta_df %>%
	filter(r_mouse_id == c_mouse_id,
		r_time_point == 'Initial' & c_time_point == 'TOC') %>%
	mutate(comparison = 'i_t')
TOC_intra_TOC <- meta_beta_df %>%
	filter(r_mouse_id != c_mouse_id,
		c_abx == r_abx,
		c_time_point == 'TOC' & r_time_point == 'TOC') %>%
	mutate(comparison = 't_t')
TOC_inter_TOC <- meta_beta_df %>%
	filter(r_mouse_id != c_mouse_id,
		c_abx != r_abx,
		c_time_point == 'TOC' & r_time_point == 'TOC') %>%
	mutate(comparison = 'tXt')
end_toC <- meta_beta_df %>%
	filter(r_mouse_id == c_mouse_id,
		r_time_point == 'End' & c_time_point == 'TOC') %>%
	mutate(comparison = 'e_t')
end_intial <- meta_beta_df %>%
	filter(r_mouse_id == c_mouse_id,
		r_time_point == 'Initial' & c_time_point == 'End') %>%
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

beta_div_df <- bind_rows(list(initial_distances, initial_inter_intial, TOC_intra_initial, TOC_intra_TOC,
	TOC_inter_TOC, end_toC, end_intial, end_intra_end, end_inter_end)) %>%
	mutate(comparison = factor(comparison,
		levels = c('i_i', 'i_t', 'e_i', 'e_t', 'iXi', 't_t', 'tXt', 'e_e', 'eXe'),
		labels = c('Initial\nvs\nInitial',
			'TOC\nvs\nInitial',
			'End\nvs\nInitial',
			'End\nvs\nTOC',
			'Initial\nvs\ninter\nInitial',
			'TOC\nvs\nintra\nTOC',
			'TOC\nvs\ninter\nTOC',
			'End\nvs\nintra\nEnd',
			'End\nvs\ninter\nEnd'))) %>%
	filter(comparison %in% c('Initial\nvs\nInitial',
		'TOC\nvs\nInitial',
		'End\nvs\nInitial',
		'End\nvs\nintra\nEnd',
		'End\nvs\ninter\nEnd'))

beta_wilcox_df <- beta_div_df %>%
	mutate(subset = paste(comparison, c_clearance, sep = '.')) %>%
	group_by(c_abx) %>%
	nest() %>%
	mutate(data = map(data, function(data) nest(group_by(data, subset)))) %>%
	mutate(data = map(data, function(nested_df){
		pvalue_df <- sapply(nested_df$data, function(x) sapply(nested_df$data, function(y)
				wilcox.test(x$distances,y$distances)$p.value)) %>%
			data.frame
		pvalue_df[upper.tri(pvalue_df)] <- NA # set all of upper triangle to 1 to eliminate duplicate comparisons
		colnames(pvalue_df) <- nested_df$subset
		pvalue_df <- pvalue_df %>%
			mutate(row_names = nested_df$subset) %>%
			pivot_longer(col = -row_names, names_to = 'col_names', values_to = 'pvalue') %>%
			filter(pvalue != 1, !is.na(pvalue)) %>% # eliminate all self comparisons and upper triangle
			mutate(pvalue = p.adjust(pvalue, method = 'BH')) # correct p values
		return(pvalue_df)
		})) %>%
	unnest(data) %>%
	ungroup

beta_signif_label_position <- data.frame(
	c_abx = c(rep('Clindamycin', 5), rep(rep(c('Cefoperazone', 'Streptomycin'), each = 5), 2)),
	subset = c(rep(paste(c('Initial\nvs\nInitial', 'TOC\nvs\nInitial', 'End\nvs\nInitial',
							'End\nvs\nintra\nEnd', 'End\nvs\ninter\nEnd'), 'Cleared', sep = '.'), 3),
			rep(paste(c('Initial\nvs\nInitial', 'TOC\nvs\nInitial', 'End\nvs\nInitial',
						'End\nvs\nintra\nEnd', 'End\nvs\ninter\nEnd'), 'Colonized', sep = '.'), 2)),
	x = c(1,2,3,1,2, rep(c(0.85, 1.85, 2.85, 0.85, 1.85), 2), rep(c(1.15, 2.15, 3.15, 1.15, 2.15), 2)),
	stringsAsFactors = F)

beta_signif_label_df <- beta_wilcox_df %>%
	filter(pvalue < 0.05) %>%
	left_join(beta_signif_label_position, by = c('c_abx', 'row_names' = 'subset')) %>%
	left_join(beta_signif_label_position, by = c('c_abx', 'col_names' = 'subset'),
		suffix = c('1','2')) %>%
	mutate(c_abx = factor(c_abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin')),
		label = '*',
		xnote = (x1 + x2)/2,
		clearance = 'NA', time_point = 'NA',
		alpha = case_when(x1 == 1.15 & x2 == 1.85 ~ 0.4,
			x2 == 1.15 & x1 == 1.85 ~ 0.4,
			x1 == 2.15 & x2 == 2.85 ~ 0.4,
			x2 == 2.15 & x1 == 2.85 ~ 0.4,
			x2 == 0.85 & x1 == 3.15 ~ 0.4,
			x1 == 0.85 & x2 == 3.15 ~ 0.4,
			x2 == 0.85 & x1 == 2.15 ~ 0.4,
			x1 == 0.85 & x2 == 2.15 ~ 0.4,
			T ~ 1))

change_from_initial <- c('Initial\nvs\nInitial', 'TOC\nvs\nInitial', 'End\nvs\nInitial')

beta_sig_initial_df <- beta_signif_label_df %>%
	filter(grepl(paste(change_from_initial, collapse = '|'), row_names),
			grepl(paste(change_from_initial, collapse = '|'), col_names)) %>%
	 mutate(y1 = c(c(1.1,1.25,1.175,1.325,1.275,1.1), # cef
	 	c(1,1.225,1.075,1.25,1.3,1.15,1), # strep
	 	c(1.1, 1.175, 1.25)), # clinda
	 	ynote = 0.01 + y1) %>%
	 filter(alpha == 1)


beta_plot <- beta_div_df %>%
	filter(comparison %in% change_from_initial) %>%
	mutate(c_abx = factor(c_abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>%
	ggplot(aes(x = comparison, y = distances, color = c_abx)) +
		geom_rect(xmin = 0, xmax = 4, ymin = 1.0000001, ymax = 2, color = NA, fill = 'white') +
		geom_point(aes(shape = c_clearance),
			position = position_jitterdodge()) +
		scale_shape_manual(values = c(1, 16), limits = c('Cleared', 'Colonized')) +
		facet_wrap(.~c_abx) +
		theme_bw() +
		labs(x = NULL, y = expression(theta[YC]), shape = 'Outcome', color = 'Time Point') +
		scale_color_manual(values = abx_color$color,
			limits = abx_color$abx) +
		guides(color = 'none') +
		scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) +
		theme(legend.position = 'bottom',
			legend.key.size = unit(0.2, 'in'),
			legend.background = element_rect(color = "black"),
			panel.spacing = unit(c(3,3),'lines'),
			axis.title = element_text(size = 14),
			strip.background = element_blank(),
			strip.text = element_blank()) +
		geom_text(data = beta_sig_initial_df,
			aes(x = xnote, y = ynote, label = label),
			show.legend = F, color = 'black') +
		geom_segment(data = beta_sig_initial_df,
			aes(x = x1, xend = x2, y = y1, yend = y1),
			size = 0.25, show.legend = F, color = 'black')

ggsave('submission/figure_2.tiff',
	plot_grid(
		plot_grid(alpha_sobs_plot + theme(), alpha_invsimpson_plot, beta_plot, ncol = 1),
		ncol = 1),
	width = 10, height = 11, units = 'in', compression = 'lzw')

ggsave('submission/figure_S3.tiff',
         plot_grid(alpha_sobs_plot_by_dose, plot_grid(NULL, alpha_invsimpson_plot_by_dose, rel_widths = c(1,40)), 
                   ncol = 1, rel_heights = c(4,5)),
       width = 8, height = 8, units = 'in', compression = 'lzw')
