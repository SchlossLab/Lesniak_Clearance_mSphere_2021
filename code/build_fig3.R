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

n_seqs <- median(apply(shared_df[,-1], 1, sum))
min_rel_abund <- 100 * 1/n_seqs

# create dataframe with difference between pre/post-colonization
delta_abund_df <- meta_df %>%
	filter(cdiff == T, clearance != 'Colonized') %>% 
	filter(day == 0 | day == last_sample) %>%  
	inner_join(shared_df, by = c('group' = 'Group')) %>% 
	pivot_longer(names_to = 'OTU', values_to = 'abundance', starts_with('Otu')) %>% 
	select(abx, mouse_id, OTU, time_point, abundance)%>% 
	group_by(mouse_id, time_point) %>% 
	mutate(total = sum(abundance),
		abundance = abundance/total * 100) %>% 
	pivot_wider(names_from = time_point, values_from = abundance)  %>% 
	rename(Beginning = 'Day 0')

test_otus <- function(antibiotic){
	tmp_shared_clearance <- delta_abund_df %>% 
		filter(abx %in% antibiotic)
	# select otus with median relative abundance > 1%
	tmp_otus <- tmp_shared_clearance %>% 
		group_by(OTU) %>% 
		gather(timepoint, abundance, Beginning, End) %>% 
		summarise(median_abundance = median(abundance, na.rm = T)) %>% 
		filter(median_abundance > 0.5) %>% 
		pull(OTU)

    warn_orig <- options("warn")$warn
	options("warn"= -1)

	sig_otus <- tmp_shared_clearance %>% 
		filter(OTU %in% tmp_otus) %>% 
		group_by(OTU) %>% 
		summarise(pvalue = wilcox.test(Beginning, End)$p.value) %>% 
		mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% 
		filter(pvalue < 0.05)

	options('warn' = warn_orig)

	tmp_shared_clearance <- tmp_shared_clearance %>%  
		inner_join(sig_otus, by = c('OTU'))  %>% 
		inner_join(tax_df, by = c('OTU')) 
	if(nrow(tmp_shared_clearance) < 1) { 
		print(paste('No significant OTUs for', paste(antibiotic, collapse = '_')))
		return(data.frame(antibiotic = paste(antibiotic, collapse = '_'), 
		stringsAsFactors = F)) 	
	} else {
	return(data.frame(tmp_shared_clearance, antibiotic = paste(antibiotic, collapse = '_'), 
		stringsAsFactors = F)) 
	}
}

treatment_list <- list('Cefoperazone', 'Clindamycin', 'Streptomycin', c('Cefoperazone', 'Clindamycin', 'Streptomycin'))

test_output_df <- map_dfr(treatment_list, ~ test_otus(.x))

differential_abundance_df <- function(df, treatment_group){
	temp_df <- df %>% 
		filter(antibiotic == paste(treatment_group, collapse = '_')) %>% 
		mutate(tax_pval = paste0(tax_otu_label, '\n(p=', 
				formatC(pvalue, format = "e", digits = 2), ')')) 
	output_df <- temp_df %>% 
		gather(time_point, abundance, Beginning, End) %>% 
		right_join(., group_by(., OTU) %>% 
				summarise(order = mean(abundance, na.rm = T)),
			by = c('OTU')) %>% 
		right_join(., group_by(., OTU, time_point) %>% 
				summarise(median_abundance = (median(abundance, na.rm = T)) + 0.04),
			by = c('OTU', 'time_point')) %>% 
		select(-abx, -mouse_id, -abundance) %>% distinct %>% 
		spread(time_point, median_abundance) %>% 
		inner_join(
			select(gather(temp_df, time_point, abundance, Beginning, End), 
				time_point, OTU, abundance))
	return(data.frame(output_df, treatment = paste(treatment_group, collapse = '_'),
		stringsAsFactors = F))
}

plot_df <- map_dfr(treatment_list, ~ differential_abundance_df(test_output_df, .x))

lod_df <- data.frame(x = 0.75, y = min_rel_abund)
diff_abund_plot <- plot_df %>% 
	filter(treatment != 'Cefoperazone_Clindamycin_Streptomycin',
		!is.na(OTU)) %>% 
	mutate(tax_pval = gsub('_unclassified', '', tax_pval),
		tax_pval = gsub('\\\n\\(p=.*\\)$', '', tax_pval)) %>% 
	ggplot(aes(x = reorder(tax_pval, -order), color = time_point)) + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 0.5, 
			linetype = 'solid', color = 'black') + 
		geom_hline(data = lod_df, aes(yintercept = y), size = 1, 
			linetype = 'dashed', color = 'white') + 
		geom_point(aes(y = (abundance) + 0.04), 
			position = position_dodge(width = 0.5), alpha = 0.2) + 
		geom_point(aes(y = Beginning), color = 'red3', size = 3) + 
		geom_point(aes(y = End), color = 'cyan4', size = 3) +
		geom_segment(aes(y = Beginning, yend = End, 
				xend = reorder(tax_pval, -order)), 
			arrow = arrow(type = 'closed', angle = 15), color = 'black') + 
		scale_y_log10(
	   		breaks = scales::trans_breaks("log10", function(x) 10^x),
	   		labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
		coord_flip() + theme_bw() +  
		labs(x = NULL, y = 'Relative Abundance (%)', color = 'Time Point') + 
		theme(legend.position = c(0.85, 0.925), 
			legend.background = element_rect(color = "black")) + 
		geom_label(data = lod_df, aes(x = x, y = y), label = "LOD", 
			fill = "white", color = 'black', label.size = NA, inherit.aes = FALSE) + 
		facet_grid(treatment~., scale = 'free_y', space = 'free') + 
		theme(text = element_text(size = 20)) + 
		guides(colour = guide_legend(override.aes = list(alpha = 1)))

ggsave('results/figures/figure_3_diff_abund_plot_long.jpg', diff_abund_plot, width = 10, height = 11, units = 'in')
#ggsave('results/figures/figure_3_diff_abund_plot_wide.jpg', diff_abund_plot, width = 15, height = 8, units = 'in')
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