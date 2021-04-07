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
	mutate(dose = factor(dose, levels = c("Clindamycin 10 mg/kg", "Cefoperazone 0.5 mg/ml",
			"Cefoperazone 0.3 mg/ml", "Cefoperazone 0.1 mg/ml", "Streptomycin 5 mg/ml", 
			"Streptomycin 0.5 mg/ml", "Streptomycin 0.1 mg/ml")))
shared_df <- read_tsv(shared_file) %>% 
	select(-label, -numOtus) %>% 
	filter(Group %in% meta_df$group)
tax_df <- read_tsv(tax_file)
source(sum_taxa_function) # function to create taxanomic labels for OTUs
	# sum_otu_by_taxa(taxonomy_df, otu_df, taxa_level = 'NA', top_n = 0, silent = T)

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))

lod_df <- data.frame(day = 2, CFU = 100, dose = '0.5 mg/ml')
# plot C difficile colonization level
plot_colonization <- function(antibiotic){
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)
	
	colonization_plot <- meta_df %>% 
		filter(day >= 0,
			abx == antibiotic,
			!is.na(CFU)) %>%
		mutate(dose = factor(gsub(paste0(antibiotic, ' '), '', dose),
				levels = c("10 mg/kg", "5 mg/ml", "0.5 mg/ml",
				"0.3 mg/ml", "0.1 mg/ml"))) %>% 
		ggplot(aes(x = day, y = CFU)) + 
			geom_line(aes(group = mouse_id), color = abx_col, alpha = 0.3) + 
			stat_summary(fun=median, geom="line", size = 1, color = abx_col) + # create median line
			scale_x_continuous(breaks = 0:10, labels = c('', '1', '', '3', '', '5', '', '7', '', '9', '')) + # make ticks for each day
			# create a dotted line labeled LOD for limit of detection
			scale_y_log10(breaks = c(10^2, 10^4, 10^6, 10^8),
					labels = c('10^2', '10^4', '10^6', '10^8')) + # scale y axis log10 and label 10^x
			theme_bw() + labs(x = 'Day', y = expression(italic('C. difficile')~' CFU')) + 
			theme(panel.grid.minor = element_blank(),
				axis.text = element_text(size = 12),
				axis.title = element_text(),
				axis.text.y = element_markdown(size = 12),
				legend.position = 'none',
				strip.background =element_rect(fill=abx_col),
				strip.text = element_text(colour = 'white', size = 12),
				plot.title = element_text(hjust = 0.5, size = 12)) + 
			facet_wrap(.~dose, nrow = 1) + 
			labs(title = antibiotic)
	if(antibiotic == 'Cefoperazone'){
		colonization_plot <- colonization_plot +
			geom_hline(yintercept = 101, linetype = 'dashed', size = 0.25) + 
			geom_label(data = lod_df, label = "LOD", fill = 'white', color = 'white') + 
			geom_text(data = lod_df, label = "LOD", color = 'black', size = 12/.pt)
	} else {
		colonization_plot <- colonization_plot +
			geom_hline(yintercept = 101, linetype = 'dashed', size = 0.25)
	}
}

# plot relative abundance over time in a heatmap plot
# mice along the x axis, taxonomic classification along the y axis and color intensity by log10 relative abundance

shared_meta_toi_df <- meta_df %>% 
	filter(day == 0)
meta_initial_df <- meta_df %>% 
	filter(time_point == 'Initial')	 %>% 
	pull(group)
toi_shared <- shared_df %>% 
	filter(Group %in% shared_meta_toi_df$group)
intial_shared <- shared_df %>% 
	filter(Group %in% meta_initial_df)

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
		facet_grid(.~dose, scales = 'free_x', space = 'free') +
		labs(x = NULL, y = NULL, #title = paste('Day 0 Community - Top', n_taxa, 'Genus'),
			fill = 'Relative Abundance (%)') + 
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
			axis.ticks.x=element_blank(), 
			strip.text = element_blank(),
			strip.background =element_blank(),
			plot.margin = margin(0, 0.1, 0, 0.1, "cm"), 
			legend.position = 'bottom',
			axis.text = element_text(size = 12)) + 
		guides(fill = guide_colorbar(title.position = "top"))
	
	outcome_plot <- abundance_df %>% 
		filter(abx == antibiotic,
			taxa == '*Clostridium sensu stricto*') %>% 
		mutate(dose = factor(gsub(paste0(antibiotic, ' '), '', dose),
				levels = c("10 mg/kg", "5 mg/ml", "0.5 mg/ml",
			"0.3 mg/ml", "0.1 mg/ml")),
			clearance = ifelse(clearance == 'Clearing', 'Colonized', clearance),
			clearance = ifelse(clearance == 'Uncolonized', 'Not Colonized', clearance),
			taxa = '     Colonization Outcome') %>% 
		ggplot(aes(x = Group, y =taxa, fill = clearance)) + 
			geom_tile() +
			theme_bw() + 
			facet_grid(.~dose, scales = 'free_x', space = 'free') +
			labs(x = NULL, y = NULL, title = NULL, fill = NULL) + 
			scale_fill_manual(limits = c('Not Colonized', 'Cleared', 'Colonized'),
				values = c('white', 'grey80', 'grey20')) +
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
				axis.ticks.x=element_blank(), 
				axis.ticks.y=element_line(colour = 'white'), 
				panel.grid=element_blank(), panel.border=element_blank(),
				legend.key=element_rect(color='black'),
				legend.key.size = unit(.5, 'line'),
				legend.text = element_text(size = 12),
				plot.margin = margin(0, 0.1, 0, 0.1, "cm"), 
				strip.background =element_rect(fill=abx_col),
				strip.text = element_text(colour = 'white', size = 12),
				plot.title = element_text(hjust = 0.5, size = 12))
	if(antibiotic == 'Cefoperazone'){
		outcome_plot <- outcome_plot +
			theme(legend.position = 'top',
				axis.text.y = element_blank())
		abundance_plot <- abundance_plot +
			theme(axis.text.y = element_blank(),
				axis.ticks.y = element_blank())
		plot_grid(outcome_plot, abundance_plot , ncol = 1, rel_heights = c(2, 10))
	} else if(antibiotic == 'Clindamycin'){
		outcome_plot <- outcome_plot +
			theme(legend.position = 'none',
				axis.text.y = element_text(size = 12))
		abundance_plot <- abundance_plot +
			theme(axis.text.y = element_markdown(size = 12))
		plot_grid(NULL, outcome_plot, abundance_plot, ncol = 1, rel_heights = c(1 , 1, 10))
	} else {
		outcome_plot <- outcome_plot +
			theme(legend.position = 'none',
				axis.text.y = element_blank())
		abundance_plot <- abundance_plot +
			theme(axis.text.y = element_blank(),
				axis.ticks.y = element_blank())
		plot_grid(NULL, outcome_plot, abundance_plot, ncol = 1, rel_heights = c(1 , 1, 10))
	}
	
}

clinda_cfu_plot <- plot_colonization('Clindamycin')
cef_cfu_plot <- plot_colonization('Cefoperazone') + 
	theme(axis.ticks.y = element_blank(),
		axis.text.y = element_blank(),
		axis.title.y = element_blank())
strep_cfu_plot <- plot_colonization('Streptomycin') + 
	theme(axis.ticks.y = element_blank(),
		axis.text.y = element_blank(),
		axis.title.y = element_blank())


clinda_abun_plot <- plot_abundance('Clindamycin')
cef_abun_plot <- plot_abundance('Cefoperazone')
strep_abun_plot <- plot_abundance('Streptomycin')
	
# plot initial community
meta_initial_df <- meta_df %>% 
	filter(time_point == 'Initial')	 %>% 
	pull(group)

intial_shared <- shared_df %>% 
	filter(Group %in% meta_initial_df)

initial_abundance_df <- sum_otu_by_taxa(tax_df, intial_shared, taxa_level = 'Genus') %>% 
	mutate(taxa = gsub('_unclassified', '', taxa),
		taxa = gsub('_', ' ', taxa),
		taxa = ifelse(taxa == 'Other', taxa, paste0('*', taxa, '*')),
		taxa = ifelse(taxa %in% taxa_order, taxa, 'Other'),
		taxa = factor(taxa, labels = taxa_order, levels = taxa_order)) %>% 
	group_by(Group, taxa) %>% 
	mutate(abundance = sum(abundance)) %>% 
	group_by(Group) %>% 
	mutate(total = sum(abundance),
		relative_abundance = log10((abundance/total * 100) + 0.01),
		day = 'Initial')

initial_abundance_plot <- initial_abundance_df %>% 
	ggplot(aes(x = Group, y =taxa, fill = relative_abundance)) + 
		geom_tile(height = 0.8) +
		scale_fill_gradient2(low="white", mid='#0B775E', high = 'black',
			limits = c(-2,2), na.value = NA, midpoint = .3,
			breaks = c(-2.5, -1, 0, 1, 2), labels = c('', '0.1', '1', '10', '100')) + 
		theme_bw() + 
		labs(x = NULL, y = NULL, #title = paste('Day 0 Community - Top', n_taxa, 'Genus'),
			fill = 'Relative Abundance (%)') + 
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
			axis.ticks.x=element_blank(), 
			axis.text.y = element_markdown(),
			panel.grid = element_blank(),
			strip.text = element_blank(),
			strip.background =element_blank(),
			plot.margin = margin(0, 0.1, 0, 0.1, "cm"), 
			legend.position = 'bottom') + 
		guides(fill = guide_colorbar(title.position = "top"))


# save plot, top row is colonization plot, middle row are diversity plots, bottom row is temporal abundance plot
ggsave('submission/figure_1.tiff', plot_grid(
		plot_grid(NULL, clinda_cfu_plot, NULL, cef_cfu_plot, NULL, strep_cfu_plot, rel_widths = c(1.5, 2.5, .25, 6, .25, 6), 
			labels = c('', 'A', '', 'B', '', 'C'), nrow = 1),
		NULL,
		plot_grid(clinda_abun_plot, NULL,  cef_abun_plot, NULL, strep_abun_plot, rel_widths = c(4, .25, 6, .25, 6), 
			labels = c('D', '', 'E', '', 'F'), nrow = 1),
		rel_heights = c(4, .1, 6), ncol = 1), 
	width = 6.87*2, height = 4.58*2, units = 'in', compression = 'lzw')

ggsave(paste0("submission/figure_S1.tiff"), 
	plot = initial_abundance_plot,
	width = 8, height = 6, units="in", compression = 'lzw')