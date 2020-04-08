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
plot_colonization_abundance <- function(antibiotic, n_taxa){
	dosages <- meta_df %>% 
		filter(abx == antibiotic) %>% 
		pull(dose) %>% 
		unique

	abx_col <- abx_color %>% 
		filter(abx == antibiotic) %>% 
		pull(color)
	abx_group <- meta_df %>% 
		filter(abx == antibiotic,
			cdiff == T,
			day == 0) %>% 
		pull(group)
	abx_shared <- shared_df %>% 
		filter(Group %in% abx_group)
	abx_meta <- meta_df %>% 
		filter(cdiff == T, day >= -1, abx == antibiotic) %>% 
		mutate(CFU = case_when(CFU == 0 ~ 60,
			T ~ CFU),
			dose = paste(antibiotic, '-', dose, 'mg/mL'),
			dose = factor(dose, levels = c(paste(antibiotic, '-', max(dosages), 'mg/mL'), 
				paste(antibiotic, '-', median(dosages), 'mg/mL') ,paste(antibiotic, '-', min(dosages), 'mg/mL')))) %>% 
		filter(!is.na(CFU))

	cfu_plot <- abx_meta %>% 
		ggplot(aes(x = day, y = CFU)) + 
			stat_summary(fun.y=median, geom="line", size = 1, color = abx_col) +
	        stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), color = abx_col) + 
	        geom_line(aes(group = mouse_id), alpha = 0.3, color = abx_col) + 
	        scale_x_continuous(breaks = -1:10) +
			annotate(x = -1, y = 200, geom = 'label', label = "LOD", 
				fill = "white", color = 'black', label.size = NA) + 
			geom_hline(yintercept = 101, linetype = 'dashed', size = 0.25) + 
			scale_y_log10(
	   			breaks = scales::trans_breaks("log10", function(x) 10^x),
	   			labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
			theme_bw() + 
			facet_wrap(.~dose, nrow = 1) + 
			labs(x = 'Day', y = expression(italic('C. difficile')~' CFU'))

	abundance_plot <- sum_otu_by_taxa(tax_df, abx_shared, taxa_level = 'Genus', top_n = n_taxa) %>% 
		left_join(select(abx_meta, group, dose), by = c('Group' = 'group')) %>% 
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
				fill = 'Day 0 Relative Abundance (%)\nColor Intesity based on Log10') + 
			theme(axis.title.x=element_blank(),
	        	axis.text.x=element_blank(),
	        	axis.ticks.x=element_blank(),
	        	axis.text.y = element_text(angle = 45),
	        	legend.position = 'bottom')
	return(plot_grid(
			plot_grid(NULL, cfu_plot, nrow = 1, rel_widths = c(1,11)),
		abundance_plot, ncol = 1))
}

cef_abundance <- plot_colonization_abundance('Cefoperazone', 12)
strep_abundance <- plot_colonization_abundance('Streptomycin', 12)

ggsave('results/figures/figure_2.jpg', plot_grid(cef_abundance, strep_abundance, labels = c('A', 'B'), ncol =1), width = 14, height = 12)