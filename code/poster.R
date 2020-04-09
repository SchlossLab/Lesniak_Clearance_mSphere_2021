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

# create a dataframe with the change in relative abundance by day
diff_rel_abund <- meta_file %>% 
	select(group, mouse_id, day) %>% 
	inner_join(rel_abund, by = 'group') %>% 
	group_by(mouse_id) %>% 
	arrange(mouse_id, day) %>% 
	mutate_at(vars(contains('Otu')), 
		function(otu) { otu - lag(otu) } ) %>% 
	ungroup

####### Colonization dynamic
cfu_lod_df <- data.frame(x = -0.5, y = 2) 
cfu_plot <- meta_file %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'), cdiff == T, day >= 0) %>% 
	ggplot(aes(x = day, y = log10CFU, group = mouse_id, color = abx)) + 
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
# input_dataframe_name <- 'cef_0.3_clearance_df' # for testing
nmds_initial_end_df <- meta_file %>% 
	filter(abx %in% c('clinda', 'cef', 'strep'), cdiff == T) %>% 
	group_by(mouse_id) %>% 
	filter(day == min_day | day == max_day)
mothur <- '/mothur/mothur'

plot_nmds <- function(input_dataframe_name){
	i <- input_dataframe_name
	current_df <- get(i)
	current_shared <- shared_file[row.names(shared_file) %in% current_df$group, ]
	# remove otus that either have 0 or only present in 2 or fewer samples
	present_otus <- current_shared %>% 
	  map_dbl(~ sum(. > 0)) %>% 
	  which(x = (. > 2)) %>% 
	  names
	current_shared <- current_shared %>% 
		mutate(label = 0.03, Group = row.names(.), 
			numOtus = length(present_otus)) %>% 
		select(label, Group, numOtus, one_of(present_otus))

	# write files to be used in mothur for nmds analysis
	write_tsv(path = paste0('data/mothur/', i, '.shared'), 
		x = current_shared)
	# run nmds
	system(paste0(mothur, ' "#set.dir(input=data/mothur, output=data/mothur);
		dist.shared(shared=', i, '.shared, calc=thetayc-jclass, subsample=t);
		nmds(phylip=', i, '.thetayc.0.03.lt.ave.dist)"'))

	# plot  results
	nmds_df <- read_tsv(paste0('data/mothur/', i, '.thetayc.0.03.lt.ave.nmds.axes'))

	nmds_plot1 <- meta_file %>% 
		inner_join(clearance, by ='mouse_id') %>% 
		filter(abx %in% c('clinda', 'strep', 'cef'), cdiff == T, clearance != 'Clearing') %>% 
		filter(day == inital_sample | day == 0 | day == last_sample) %>% 
		mutate(pt_size = ifelse(time_point == "End", 3, 1),
			pt_shape = ifelse(time_point == 'End', 3, 1),
			clearance = case_when(clearance == 'Cleared' ~'Cleared Colonization',
				clearance == 'Colonized' ~ 'Remain Colonized')) %>% 
		inner_join(nmds_df) %>% 
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
	ggsave('~/Desktop/nmds_temporal_initial_end_by_abx.jpg', nmds_plot,
			width = 6, height = 4)
}
plot_nmds('nmds_initial_end_df')
######## Community changes
#### 
# whats the difference between the communities within antibiotic between day 0 and end point
####

delta_abund_df <- meta_file %>%
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'), cdiff == T, 
		clearance != 'Colonized'#, 
		#mouse_id %in% c('42_1', '600_2', '88_1') #mice were being excluded since they were missing either day 0 or day 10
		# but removing and adding a filtering step at the beginning of the script to only include samples in both shared/meta
		) %>% 
	filter(day == 0 | day == last_sample) %>%  
	inner_join(rel_abund, by = 'group') %>% 
	gather(OTU, abundance, contains('Otu00')) %>% 
	select(abx, mouse_id, OTU, time_point, abundance) %>% 
	spread(time_point, abundance) %>% 
	rename(start = 'Day 0', end = End)
# dont need to remove since not comparing individuals #filter( !( is.na(start) | is.na(end) ) ) %>% # remove any mice missing either sample 
# dont need to remove since not comparing individuals #filter(end == 0 & start == 0) 

test_otus <- function(antibiotic){
	tmp_shared_clearance <- delta_abund_df %>% 
		filter(abx %in% antibiotic)
	# select otus with median relative abundance > 1%
	tmp_otus <- tmp_shared_clearance %>% 
		group_by(OTU) %>% 
		gather(timepoint, abundance, end, start) %>% 
		summarise(median_abundance = median(abundance, na.rm = T)) %>% 
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
	if(nrow(tmp_shared_clearance) < 1) { 
		print(paste('No significant OTUs for', paste(antibiotic, collapse = '_')))
		return(data.frame(antibiotic = paste(antibiotic, collapse = '_'), 
		stringsAsFactors = F)) 	
	} else {
	return(data.frame(tmp_shared_clearance, antibiotic = paste(antibiotic, collapse = '_'), 
		stringsAsFactors = F)) 
	}
}

treatment_list <- list('Cefoperazone', 'Clindamycin', 'Streptomycin', 
	c('Cefoperazone', 'Clindamycin', 'Streptomycin'))

test_output_df <- map_dfr(treatment_list, ~ test_otus(.x))

differential_abundance_df <- function(df, treatment_group){
	temp_df <- df %>% 
		filter(antibiotic == paste(treatment_group, collapse = '_')) %>% 
		mutate(tax_pval = paste0(tax_otu_label, '\n(p=', 
				formatC(pvalue, format = "e", digits = 2), ')')) 
	output_df <- temp_df %>% 
		gather(time_point, abundance, start, end) %>% 
		right_join(., group_by(., OTU) %>% 
				summarise(order = mean(abundance, na.rm = T)),
			by = c('OTU')) %>% 
		right_join(., group_by(., OTU, time_point) %>% 
				summarise(median_abundance = (median(abundance, na.rm = T)) + 0.04),
			by = c('OTU', 'time_point')) %>% 
		select(-abx, -mouse_id, -abundance) %>% distinct %>% 
		spread(time_point, median_abundance) %>% 
		inner_join(
			select(gather(temp_df, time_point, abundance, start, end), 
				time_point, OTU, abundance))
	return(data.frame(output_df, treatment = paste(treatment_group, collapse = '_'),
		stringsAsFactors = F))
}

plot_df <- map_dfr(treatment_list, ~ differential_abundance_df(test_output_df, .x))

lod_df <- data.frame(x = 0.75, y = min_rel_abund)
diff_abund_plot <- plot_df %>% 
	filter(treatment != 'cef_clinda_strep',
		!is.na(OTU)) %>% 
	mutate(tax_pval = gsub('_unclassified', '', tax_pval)) %>% 
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