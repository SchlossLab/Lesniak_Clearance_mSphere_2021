##############
#
# run script to generate plots for Figure 2
#	Can other conditions naturally clear C. difficile colonization?
# 
# Nick Lesniak 04-06-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/sample.final.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#
##############



library(tidyverse)
library(cowplot)
library(ggtext)  # remotes::install_github("wilkelab/ggtext")

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
alpha_div_file <- 'data/mothur/sample.final.groups.ave-std.summary'
beta_div_file <- 'data/mothur/sample.final.thetayc.0.03.lt.ave.dist'
dist_function <- 'code/read_dist.R'

# read in data
meta_df   <- read_tsv(meta_file) %>% 
	filter(clearance %in% c('Cleared', 'Colonized'),
		cdiff == T, time_point != 'Intermediate') %>%  
	mutate(time_point = ifelse(time_point == 'Day 0', 'TOI', time_point),
		time_point = factor(time_point, levels = c('Initial', 'TOI', 'End'),
			labels = c('Initial', 'TOI', 'End')),
		abx = factor(abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'),
			labels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>% 
	select(group, dose, day, abx, mouse_id, clearance, time_point) 

alpha_df <- read_tsv(alpha_div_file) %>% 
	filter(group %in% meta_df$group,
		method == 'ave')

source(dist_function) # function to read in distance file and convert from triangle to dataframe
beta_df <- read_dist(beta_div_file) %>% 
	filter(rows %in% meta_df$group,
		columns %in% meta_df$group)

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))

# plot Sobs by day
alpha_sobs_plot <- alpha_df %>% 
	select(group, sobs) %>% 
	inner_join(select(meta_df, abx, group, time_point, clearance), by = c('group')) %>% 
	ggplot(aes(x = time_point, y = sobs, color = abx, fill = clearance, group = interaction(clearance, time_point))) + 
		geom_boxplot(position = 'dodge') + 
		#geom_point(alpha = 0.4, position = position_jitterdodge()) + 
		coord_cartesian(ylim = c(0,160)) +
        #scale_x_continuous(breaks = -1:10) +
		scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) +
		scale_color_manual(values = abx_color$color, limits = abx_color$abx) + 
		theme_bw() + labs(x = 'Day', y = expression(~S[obs])) + 
		theme(panel.grid.minor = element_blank(),
			legend.position = 'none',
			panel.spacing = unit(c(3,3),'lines')) + 
		facet_wrap(.~abx)
# add labels to plots for figure
alpha_sobs_plot <- plot_grid(
	plot_grid(NULL, NULL, NULL, labels = c('A', 'B', 'C'), nrow = 1), 
	alpha_sobs_plot, ncol = 1, rel_heights = c(1,19))


# plot inverse simpson by day
alpha_invsimpson_plot <- alpha_df %>% 
	select(group, invsimpson) %>% 
	inner_join(select(meta_df, abx, group, time_point, clearance), by = c('group')) %>% 
	ggplot(aes(x = time_point, y = invsimpson, color = abx, fill = clearance, group = interaction(clearance, time_point))) + 
		geom_boxplot(position = 'dodge') + 
		#geom_point(alpha = 0.4, position = position_jitterdodge()) + 
		#coord_cartesian(ylim = c(0,160)) +
        #scale_x_continuous(breaks = -1:10) +
		scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) +
		scale_color_manual(values = abx_color$color, limits = abx_color$abx) + 
		theme_bw() + labs(x = 'Day', y = 'Inverse Simpson') + 
		theme(panel.grid.minor = element_blank(),
			legend.position = 'none',
			panel.spacing = unit(c(3,3),'lines')) + 
		facet_wrap(.~abx)

###############################################################################
#   How different are communities at TOI and End?
####

# plot Theta yc differences for 
#	initial as control range, 
#	TOI vs own initial, TOI vs other Initial, 
#	End vs End, End vs Initial, End vs TOI, End vs other End
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
TOI_intra_initial <- meta_beta_df %>% 
	filter(r_mouse_id == c_mouse_id,
		r_time_point == 'Initial' & c_time_point == 'TOI') %>% 
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
		r_time_point == 'End' & c_time_point == 'TOI') %>% 
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

beta_div_df <- bind_rows(list(initial_distances, initial_inter_intial, TOI_intra_initial, TOI_intra_TOI, 
	TOI_inter_TOI, end_toi, end_intial, end_intra_end, end_inter_end)) %>% 
	mutate(comparison = factor(comparison, 
		levels = c('i_i', 'i_t', 'e_i', 'e_t', 'iXi', 't_t', 'tXt', 'e_e', 'eXe'),
		labels = c('Initial\nvs\nInitial', 
			'TOI\nvs\nInitial', 
			'End\nvs\nInitial', 
			'End\nvs\nTOI', 
			'Initial\nvs\ninter\nInitial',
			'TOI\nvs\nintra\nTOI', 
			'TOI\nvs\ninter\nTOI',
			'End\nvs\nintra\nEnd', 
			'End\nvs\ninter\nEnd')))

beta_plot <- beta_div_df %>% 
	filter(comparison %in% c('Initial\nvs\nInitial', 
			'TOI\nvs\nInitial', 
			'End\nvs\nInitial')) %>% 
	mutate(c_abx = factor(c_abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>% 
	ggplot(aes(x = comparison, y = distances, color = c_abx)) + 
		coord_cartesian(ylim = c(0,1)) +
		geom_boxplot(aes(fill = c_clearance)) + 
		scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) +
		facet_wrap(.~c_abx) + 
		theme_bw() + 
		labs(x = NULL, y = expression(theta[YC]), fill = 'Outcome', color = 'Time Point') + 
		scale_color_manual(values = abx_color$color,
			limits = abx_color$abx) + 
		theme(legend.position = 'bottom',
			legend.key.size = unit(0.2, 'in'),
			legend.background = element_rect(color = "black"),
			panel.spacing = unit(c(3,3),'lines'))

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


ggsave('results/figures/figure_2.jpg', 
	plot_grid(
		plot_grid(alpha_sobs_plot, alpha_invsimpson_plot, beta_plot, ncol = 1),
		ncol = 1), 
	width = 10, height = 15, units = 'in')

ggsave('results/figures/figure_S1.jpg', beta_supp_plot, 
	width = 10, height = 10, units = 'in')