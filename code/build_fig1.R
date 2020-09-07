##############
#
# run script to generate plots for Figure 1
#	What occurs while C. difficile colonization is naturally cleared?
# 
# Nick Lesniak 04-06-2020
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
library(ggtext)


meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
shared_file <- 'data/mothur/sample.final.0.03.subsample.shared'
sum_taxa_function <- 'code/sum_otu_by_taxa.R'

# read in data
meta_df <- read_tsv(meta_file) %>% 
	filter(cdiff == T) %>% 
	mutate(CFU = case_when(CFU == 0 ~ 60, T ~ CFU), # shift 0 counts to just below limit of detection line
		dose = ifelse(abx == 'Clindamycin', paste(abx, dose, 'mg/kg'),
			paste(abx, dose, 'mg/ml'))) %>% 
	group_by(dose) %>% 
	mutate(n = length(unique(mouse_id))) %>% 
	ungroup() %>% 
	mutate(dose = paste(dose, '(N =', n, ')'),
		dose = factor(dose, levels = c("Clindamycin 10 mg/kg (N = 11 )", "Cefoperazone 0.5 mg/ml (N = 6 )",
			"Cefoperazone 0.3 mg/ml (N = 13 )", "Cefoperazone 0.1 mg/ml (N = 6 )", "Streptomycin 5 mg/ml (N = 8 )", 
			"Streptomycin 0.5 mg/ml (N = 9 )", "Streptomycin 0.1 mg/ml (N = 11 )")))
shared_df <- read_tsv(shared_file) %>% 
	select(-label, -numOtus) %>% 
	filter(Group %in% meta_df$group)
tax_df <- read_tsv(tax_file)
source(sum_taxa_function) # function to create taxanomic labels for OTUs
	# sum_otu_by_taxa(taxonomy_df, otu_df, taxa_level = 'NA', top_n = 0, silent = T)

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))

lod_df <- data.frame(day = 2, CFU = 100, dose = 'Clindamycin 10 mg/kg (N = 11 )')
# plot C difficile colonization level
colonization_plot <- meta_df %>% 
	filter(day >= 0,
		!is.na(CFU)) %>%
	ggplot(aes(x = day, y = CFU, color = abx)) + 
		geom_line(aes(group = mouse_id), alpha = 0.3) + 
		stat_summary(fun=median, geom="line", size = 1) + # create median line
		scale_color_manual(limits = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
			values = c('#D37A1F', '#3A9CBC', '#A40019')) + 
		scale_x_continuous(breaks = -1:10) + # make ticks for each day
		# create a dotted line labeled LOD for limit of detection
		geom_hline(yintercept = 101, linetype = 'dashed', size = 0.25) + 
		geom_label(data = lod_df, label = "LOD", fill = 'white', color = 'white') + 
		geom_text(data = lod_df, label = "LOD", color = 'black') + 
		scale_y_log10(breaks = c(10^2, 10^4, 10^6, 10^8),
				labels = c('10^2', '10^4', '10^6', '10^8')) + # scale y axis log10 and label 10^x
		theme_bw() + labs(x = 'Day', y = expression(italic('C. difficile')~' CFU')) + 
		theme(panel.grid.minor = element_blank(),
			axis.text.y = element_markdown(),
			legend.position = 'none',
			panel.spacing = unit(c(8,1,1,8,1,1),'lines')) + 
		facet_wrap(.~dose, nrow = 1)


# plot relative abundance over time in a heatmap plot
# mice along the x axis, taxonomic classification along the y axis and color intensity by log10 relative abundance

shared_meta_toi_df <- meta_df %>% 
	filter(day == 0)
toi_shared <- shared_df %>% 
	filter(Group %in% shared_meta_toi_df$group)

taxa_order <- rev(c("*Enterobacteriaceae*", "*Lactobacillus*", "Other", "*Porphyromonadaceae*", 
	"*Bacteroides*", "*Akkermansia*", "*Barnesiella*", "*Lachnospiraceae*", 
	"*Ruminococcaceae*", "*Pseudomonas*","*Alistipes*", "*Clostridiales*", 
	"*Clostridium sensu stricto*"))

abundance_df <- sum_otu_by_taxa(tax_df, toi_shared, taxa_level = 'Genus', top_n = 12) %>% 
	left_join(select(shared_meta_toi_df, abx, group, dose, clearance), by = c('Group' = 'group')) %>% 
	group_by(Group) %>% 
	mutate(total = sum(abundance),
		relative_abundance = log10((abundance/total * 100) + 0.01),
		taxa = gsub('_unclassified', '', taxa),
		taxa = gsub('_', ' ', taxa),
		taxa = ifelse(taxa == 'Other', taxa, paste0('*', taxa, '*')),
		taxa = factor(taxa, labels = taxa_order, levels = taxa_order),
		day = 'Time of\nInfection')

plot_abundance <- function(antibiotic){
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)
	
	abundance_plot <- abundance_df %>% 
		filter(abx == antibiotic) %>% 
	ggplot(aes(x = Group, y =taxa, fill = relative_abundance)) + 
		geom_tile(height = 0.8) +
		scale_fill_gradient2(low="white", mid=abx_col, high = 'black',
			limits = c(-2,2), na.value = NA, midpoint = .3,
			breaks = c(-2.5, -1, 0, 1, 2), labels = c('', '0.1', '1', '10', '100')) + 
		theme_bw() + 
		facet_wrap(dose~., scales = 'free_x', nrow = 1) +
		labs(x = NULL, y = NULL, #title = paste('Day 0 Community - Top', n_taxa, 'Genus'),
			fill = 'Relative Abundance (%)') + 
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
			axis.ticks.x=element_blank(), 
			axis.text.y = element_markdown(),
			strip.text = element_blank(),
			strip.background =element_blank(),
			plot.margin = margin(0, 0.1, 0, 0.1, "cm"), 
			legend.position = 'bottom') + 
		guides(fill = guide_colorbar(title.position = "top"))
	
	outcome_plot <- abundance_df %>% 
		filter(abx == antibiotic,
			taxa == '*Clostridium sensu stricto*') %>% 
		mutate(dose = factor(gsub(paste0(antibiotic, ' '), '', dose),
				levels = c("10 mg/kg (N = 11 )", "0.5 mg/ml (N = 6 )",
			"0.3 mg/ml (N = 13 )", "0.1 mg/ml (N = 6 )", "5 mg/ml (N = 8 )", 
			"0.5 mg/ml (N = 9 )", "0.1 mg/ml (N = 11 )"))) %>% 
	ggplot(aes(x = Group, y =taxa, fill = clearance)) + 
		geom_tile() +
		theme_bw() + 
		facet_wrap(dose~., scales = 'free_x', nrow = 1) +
		labs(x = NULL, y = NULL, title = antibiotic,
			fill = 'Outcome') + 
		scale_fill_manual(limits = c('Uncolonized', 'Cleared', 'Clearing', 'Colonized'),
			values = c('white', 'grey80', 'grey50', 'grey20')) +
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
			axis.ticks.x=element_blank(), 
			axis.text.y=element_markdown(colour = 'white'),
			axis.ticks.y=element_line(colour = 'white'), 
			panel.grid=element_blank(), panel.border=element_blank(),
			legend.position = 'none',
			plot.margin = margin(0, 0.1, 0, 0.1, "cm"), 
			strip.background =element_rect(fill=abx_col),
			strip.text = element_text(colour = 'white'),
			plot.title = element_text(hjust = 0.5),
			)
	plot_grid(outcome_plot, abundance_plot, ncol = 1, rel_heights = c(2, 12))
}

clinda_abun_plot <- plot_abundance('Clindamycin')
cef_abun_plot <- plot_abundance('Cefoperazone')
strep_abun_plot <- plot_abundance('Streptomycin')


# save plot, top row is colonization plot, middle row are diversity plots, bottom row is temporal abundance plot
ggsave('results/figures/figure_1.jpg', plot_grid(
		plot_grid(
			plot_grid(NULL, NULL, NULL, labels = c('A', 'B', 'C'), rel_widths = c(2, 3, 3), nrow = 1),
			plot_grid(NULL, colonization_plot, rel_widths = c(1, 16)),
			rel_heights = c(1, 15), ncol = 1),
		plot_grid(
			plot_grid(NULL, NULL, NULL, labels = c('D', 'E', 'F'), rel_widths = c(2, 3, 3), nrow = 1),
			plot_grid(clinda_abun_plot, cef_abun_plot, strep_abun_plot, rel_widths = c(1.8, 4, 4), nrow = 1), 
			rel_heights = c(1, 15), ncol = 1),
	ncol = 1), width = 15, height = 10)
