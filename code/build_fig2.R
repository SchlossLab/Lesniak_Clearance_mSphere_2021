##############
#
# run script to generate plots for Figure 2
#	Can other conditions naturally clear C. difficile colonization?
# 
# Nick Lesniak 04-06-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#
##############


library(tidyverse)
library(cowplot)
library(ggtext)

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
sum_taxa_function <- 'code/sum_otu_by_taxa.R'

abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))

# read in data
meta_df   <- read_tsv(meta_file) %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'))
shared_df <- read_tsv(shared_file) %>% 
	select(-label, -numOtus) %>% 
	filter(Group %in% meta_df$group)
tax_df <- read_tsv(tax_file)
source(sum_taxa_function) # function to create taxanomic labels for OTUs
	# sum_otu_by_taxa(taxonomy_df, otu_df, taxa_level = 'NA', top_n = 0, silent = T){

# plot colonization dynamics and Day 0 abundances for Strep and Cef, faceted by dosage of antibiotic
plot_colonization_abundance <- function(antibiotic, n_taxa, label_input){
	# get vector of antibiotic dosages for labeling
	dosages <- meta_df %>% 
		filter(abx == antibiotic) %>% 
		pull(dose) %>% 
		unique
	# select color for antibiotic
	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)
	# create subset shared with only antibiotic and day 0
	abx_group <- meta_df %>% 
		filter(abx == antibiotic,
			cdiff == T,
			day == 0) %>% 
		pull(group)
	abx_shared <- shared_df %>% 
		filter(Group %in% abx_group)
	# filter metadata and modify dosage for plot labeling
	abx_meta <- meta_df %>% 
		filter(cdiff == T, day >= -1, abx == antibiotic) %>% 
		mutate(CFU = case_when(CFU == 0 ~ 60, T ~ CFU),
			dose = paste(antibiotic, '-', dose, 'mg/mL'),
			dose = factor(dose, levels = c(paste(antibiotic, '-', max(dosages), 'mg/mL'), 
				paste(antibiotic, '-', median(dosages), 'mg/mL') ,paste(antibiotic, '-', min(dosages), 'mg/mL')))) %>% 
		filter(!is.na(CFU))
	# plot cfu over time, with median/IQR bars overall plus lines for individual mice
	lod_label_df <- data.frame(day = 8, CFU = 90, fill = 'white', color = 'black',
		dose = factor(levels(abx_meta$dose)[1], levels = levels(abx_meta$dose)))
	cfu_plot <- abx_meta %>% 
		filter(day >= 0) %>% 
		ggplot(aes(x = day, y = CFU)) + 
			stat_summary(fun=median, geom="line", size = 1, color = abx_col) +
	        geom_line(aes(group = mouse_id), alpha = 0.3, color = abx_col) + 
	        scale_x_continuous(breaks = -1:10) +
			geom_hline(yintercept = 90, linetype = 'dashed', size = 0.25) + 
			geom_label(data = lod_label_df, label = "LOD", color = 'white') + 
			geom_text(data = lod_label_df, label = "LOD") + 
			scale_y_log10(breaks = c(10^2, 10^4, 10^6, 10^8),
				labels = c('10^2', '10^4', '10^6', '10^8')) + 
			theme_bw() + 
			facet_wrap(.~dose, nrow = 1) + 
			labs(x = 'Day', y = expression(italic('C. difficile')~' CFU')) + 
			theme(panel.grid.minor = element_blank(),
				panel.spacing = unit(c(3,1),'lines'),
				axis.text.y = element_markdown())
	# plot day 0 relative abundance by antibiotic dosage
	if(antibiotic == 'Cefoperazone'){
			taxa_order <- c('*Lactobacillus*', 'Other', 
				'*Escherichia/Shigella*', '*Clostridium sensu stricto*', 
				'*Enterobacteriaceae*', '*Bacteroidales*', 
				'*Alistipes*', '*Ruminococcaceae*', '*Barnesiella*', '*Akkermansia*', 
				'*Lachnospiraceae*', '*Porphyromonadaceae*', '*Bacteroides*')
		} else if(antibiotic == 'Streptomycin'){
			taxa_order <- c('*Porphyromonadaceae*','*Bacteroides*', 
				'*Akkermansia*', '*Barnesiella*', '*Bacteroidales*', 
				'*Lachnospiraceae*', '*Lactobacillus*', '*Ruminococcaceae*', 
				'*Olsenella*', '*Alistipes*', '*Anaeroplasma*', 
				'Other', '*Oscillibacter*')
	}

	abundance_plot <- sum_otu_by_taxa(tax_df, abx_shared, taxa_level = 'Genus', top_n = n_taxa) %>% 
		left_join(select(abx_meta, group, dose), by = c('Group' = 'group')) %>% 
		group_by(Group) %>% 
		mutate(total = sum(abundance),
			relative_abundance = abundance/total * 100,
			taxa = gsub('_unclassified', '', taxa),
			taxa = gsub('_', ' ', taxa),
			taxa = ifelse(taxa == 'Other', taxa, paste0('*', taxa, '*')),
			taxa = factor(taxa, levels = taxa_order),
			day = 'Time of\nInfection') %>% 
		group_by(dose, taxa, day) %>% 
		summarise(relative_abundance = log10(mean(relative_abundance) + 0.01)) %>% 
		ggplot(aes(x = day, y =taxa, fill = relative_abundance)) + 
			geom_tile(height = 0.8) +
			scale_fill_gradient2(low="white", mid=abx_col, high = 'black',
				limits = c(-2,1.8), na.value = NA, midpoint = .3,
				breaks = c(-2.5, -1, 0, 1, 2), labels = c('', '0.1', '1', '10', '100')) + 
			theme_bw() + 
			facet_wrap(dose~., scales = 'free_x', nrow = 1) +
			labs(x = NULL, y = NULL, #title = paste('Day 0 Community - Top', n_taxa, 'Genus'),
				fill = expression('Color Intesity based on Log'[10]*' Mean Relative Abundance (%)')) + 
			theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
	        	axis.ticks.x=element_blank(), 
	        	axis.text.y = element_markdown(),
	        	legend.position = 'bottom')
	# return a plot with top row colonization plot shifted right to align with abundance plot
	if(label_input){
			return(plot_grid(
				plot_grid(NULL, NULL, NULL, nrow = 1, rel_widths = c(1,2,3), labels = c('A', 'B', 'C')),
				plot_grid(cfu_plot, abundance_plot, nrow = 1),
				ncol = 1, rel_heights = c(1, 20)))
		}else{
			return(plot_grid(cfu_plot, abundance_plot, nrow = 1))
	}
}


# create plots with the top 12 genus for relative abundance plot
cef_abundance <- plot_colonization_abundance('Cefoperazone', 12, T)
strep_abundance <- plot_colonization_abundance('Streptomycin', 12, F)

ggsave('results/figures/figure_2.jpg', 
	plot_grid(cef_abundance, strep_abundance, ncol =1), 
	width = 14, height = 12)