library(tidyverse)
library(cowplot)

# file names relative to code directory for Rmd
meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
sum_tax_function <- 'code/sum_otu_by_taxa.R'

# read in data
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	mutate(mouse_id = paste(cage, mouse, sep = '_'))
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F, row.names = "Group") %>% 
	select(-label, -numOtus)
tax_df <- read.table(tax_file, sep = '\t', header = T, stringsAsFactors = F)
source(sum_taxa_function) # function to create taxanomic labels for OTUs

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
	mutate(delta_trend = ifelse(sign(delta_cfu) == -1, -1, 1))
# create a dataframe with the change in relative abundance by day
diff_rel_abund <- meta_file %>% 
	select(group, mouse_id, day) %>% 
	inner_join(rel_abund, by = 'group') %>% 
	group_by(mouse_id) %>% 
	arrange(mouse_id, day) %>% 
	mutate_at(vars(contains('Otu')), 
		function(otu) { otu - lag(otu) } ) %>% 
	ungroup

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
		inital_sample = min(day),
		max_cfu = max(CFU)) %>% 
	left_join(differenced_cfu) %>% 
	filter(day == last_sample, max_cfu > 0) %>% # remove mice that don't become colonized at all
	select(mouse_id, end_point_cfu = CFU, delta_trend, delta_cfu, last_sample, inital_sample) %>% 
	mutate(clearance = case_when(end_point_cfu == 0 ~ 'Cleared', 
		delta_trend < 0 & end_point_cfu < 100000 ~ 'Clearing', # 10^5 separates the mice 
		T ~ 'Colonized'))

####### Colonization dynamic
cfu_lod_df <- data.frame(x = -0.5, y = 2) 
cfu_plot <- meta_file %>% 
	filter(abx %in% c('clinda', 'strep', 'cef'), cdiff == T, day >= 0) %>% 
	mutate(abx = case_when(abx == 'clinda' ~ 'Clindamycin',
		abx == 'strep' ~ 'Streptomycin',
		abx == 'cef' ~ 'Cefoperazone')) %>% 
	ggplot(aes(x = day, y = log10(CFU +60), group = mouse_id, color = abx)) + 
		geom_line() + 
		facet_wrap(.~abx, ncol = 1) + 
		theme_bw() + theme(legend.position = 'none') + 
		labs(x = 'Day', y = 'C. difficile CFU (log10)') +
		geom_hline(yintercept = 2, linetype = 'dashed', size = 0.25) + 
		geom_label(data = cfu_lod_df, aes(x = x, y = y), label = "LOD", 
			fill = "white", color = 'black', label.size = NA, inherit.aes = FALSE) + 
		scale_x_continuous(breaks = 0:10) + scale_y_continuous(breaks = 2:8) +
		theme(panel.grid.minor = element_blank(), 
			text = element_text(size = 20))




####### Beta diversity

time_nmds_2m <- read.table('data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.nmds.2.axes', 
  sep = '\t', header = T)

nmds_plot <- meta_file %>% 
	inner_join(clearance, by ='mouse_id') %>% 
	filter(abx %in% c('clinda', 'strep', 'cef'), cdiff == T, clearance != 'Clearing') %>% 
	filter(day == inital_sample | day == 0 | day == last_sample) %>% 
	mutate(abx = case_when(abx == 'clinda' ~ 'Clindamycin',
			abx == 'strep' ~ 'Streptomycin',
			abx == 'cef' ~ 'Cefoperazone'),
		time_point = factor(case_when(day == 0 ~ 'Day 0',
			day < 0 ~ 'Initial',
			day > 0 ~ 'End'), levels = c('Initial', 'Day 0', 'End')),
		pt_size = ifelse(time_point == "End", 3, 1),
		pt_shape = ifelse(time_point == 'End', 3, 1),
		clearance = case_when(clearance == 'Cleared' ~'Cleared Colonization',
			clearance == 'Colonized' ~ 'Remain Colonized')) %>% 
	inner_join(time_nmds_2m) %>% 
	arrange(mouse_id, day) %>% 
	ggplot(aes(x = axis1, y = axis2, color = abx)) +
		geom_point(color = 'black', size = 0.5) + 
		geom_path(aes(group = mouse_id), 
			  arrow = arrow(type = "closed", angle = 15, length = unit(0.25, "inches"))) + 
		theme_bw() +
		labs(x = 'NMDS Axis 1', y = 'NMDS Axis 2') + 
		facet_grid(.~clearance) + 
		theme(legend.position = 'none', panel.grid.minor = element_blank(), 
			text = element_text(size = 20))


######## Community changes
#### 
# whats the difference between the communities within antibiotic between day 0 and end point
####

delta_abund_df <- meta_file %>% 
	inner_join(clearance, by = 'mouse_id') %>% 
	filter(abx %in% c('clinda', 'strep', 'cef'), cdiff == T, 
		clearance != 'Colonized', !mouse_id %in% c('42_1', '600_2', '88_1') ) %>% 
	filter(day == 0 | day == last_sample) %>% 
	inner_join(rel_abund, by = 'group') %>% 
	gather(OTU, abundance, contains('Otu00')) %>% 
	mutate(day = ifelse(day == 0, 'start', 'end')) %>% 
	select(abx, mouse_id, OTU, day, abundance) %>% 
	spread(day, abundance) %>% 
	group_by(mouse_id, OTU) %>% 
	mutate(delta_abund = end - start) %>% 
	filter(delta_abund != 0)

test_otus <- function(antibiotic, time_point){
	tmp_shared_clearance <- delta_abund_df %>% 
		filter(abx %in% antibiotic)

	# select otus with median relative abundance > 1%
	tmp_otus <- tmp_shared_clearance %>% 
		group_by(OTU) %>% 
		gather(timepoint, abundance, end, start) %>% 
		summarise(median_abundance = median(abundance)) %>% 
		filter(median_abundance > 0.5) %>% 
		pull(OTU)

    warn_orig <- options("warn")$warn
	options("warn"= -1)

	sig_otus <- tmp_shared_clearance %>% 
		filter(OTU %in% tmp_otus) %>% 
		group_by(OTU) %>% 
		summarise(pvalue = wilcox.test(end, start)$p.value) %>% 
		mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% 
		filter(pvalue < 0.05)

	options('warn' = warn_orig)

	tmp_shared_clearance <- tmp_shared_clearance %>%  
		inner_join(sig_otus, by = c('OTU'))  %>% 
		inner_join(tax_df, by = c('OTU')) 

	return(data.frame(tmp_shared_clearance, antibiotic = paste(antibiotic, collapse = '_'), 
		stringsAsFactors = F)) 
}

treatment_list <- list('cef', 'clinda', 'strep', c('cef', 'clinda', 'strep'))

test_output_df <- map_dfr(treatment_list, ~ test_otus(.x))

differential_abundance_df <- function(df, treatment_group){
	lod_df <- data.frame(x = 0.75, y = min_rel_abund)
	temp_df <- df %>% 
		filter(antibiotic == paste(treatment_group, collapse = '_')) %>% 
		mutate(tax_pval = paste0(tax_otu_label, '\n(p=', 
				formatC(pvalue, format = "e", digits = 2), ')')) 
	output_df <- temp_df %>% 
		gather(time_point, abundance, start, end) %>% 
		right_join(., group_by(., OTU) %>% 
				summarise(order = mean(abundance)),
			by = c('OTU')) %>% 
		right_join(., group_by(., OTU, time_point) %>% 
				summarise(median_abundance = (median(abundance)) + 0.04),
			by = c('OTU', 'time_point')) %>% 
		select(-abx, -mouse_id, -delta_abund, -abundance) %>% distinct %>% 
		spread(time_point, median_abundance) %>% 
		inner_join(
			select(gather(temp_df, time_point, abundance, start, end), 
				time_point, OTU, abundance))
	return(data.frame(output_df, treatment = paste(treatment_group, collapse = '_'),
		stringsAsFactors = F))
}

plot_df <- map_dfr(treatment_list, ~ differential_abundance_df(test_output_df, .x))

diff_abund_plot <- plot_df %>% 
	filter(treatment != 'cef_clinda_strep') %>% 
	mutate(treatment = case_when(treatment == 'clinda' ~ 'Clindamycin',
		treatment == 'strep' ~ 'Streptomycin',
		treatment == 'cef' ~ 'Cefoperazone'),
	tax_pval = gsub('_unclassified', '', tax_pval),
	time_point = case_when(time_point == 'start' ~ 'Day 0',
		time_point == 'end' ~ 'Day 10')) %>% 
	ggplot(aes(x = reorder(tax_pval, -order), color = time_point)) + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
			linetype = 'solid', color = 'black') + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
			linetype = 'dashed', color = 'white') + 
		geom_point(aes(y = (abundance) + 0.04), 
			position = position_dodge(width = 0.5), alpha = 0.2) + 
		geom_point(aes(y = start), color = 'red3', size = 3) + 
		geom_point(aes(y = end), color = 'cyan4', size = 3) +
		geom_segment(aes(y = start, yend = end, 
				xend = reorder(tax_pval, -order)), 
			arrow = arrow(type = 'closed', angle = 15), color = 'black') + 
		scale_y_log10() + coord_flip() + 
		theme_bw() +  
		labs(x = NULL, y = 'Relative Abundance (%)', color = 'Time Point') + 
		theme(legend.position = c(0.85, 0.925), 
			legend.background = element_rect(color = "black")) + 
		geom_label(data = lod_df, aes(x = x, y = y), label = "LOD", 
			fill = "white", color = 'black', label.size = NA, inherit.aes = FALSE) + 
		facet_grid(treatment~., scale = 'free_y', space = 'free') + 
		theme(text = element_text(size = 20)) + 
		guides(colour = guide_legend(override.aes = list(alpha = 1)))

ggsave('~/Desktop/diff_abund_plot.pdf', diff_abund_plot, width = 10, height = 11, units = 'in')
ggsave('~/Desktop/cfu_plot.pdf', cfu_plot, width = 9, height = 9, units = 'in')
ggsave('~/Desktop/nmds_plot.pdf', nmds_plot, width = 9, height = 9, units = 'in')