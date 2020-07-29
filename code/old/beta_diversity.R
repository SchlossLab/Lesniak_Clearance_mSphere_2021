library(tidyverse)
library(cowplot)

# file names relative to code directory for Rmd
meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
dist_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist'
dist_function <- 'code/read.dist.R'

# read in data
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone')) %>% 
	mutate(time_point = factor(time_point, 
			levels = c('Initial', 'Day 0', 'Intermediate', 'End')),
		dose_level = factor(dose_level, 
			levels = c('Low', 'Mid', 'High')))
source(dist_function)
tyc_df <- read_dist(dist_file)

# for control variation, determine variation between samples on initial day
initial_intra_cage <- meta_file %>% 
	filter(time_point == 'Initial') %>% 
	select(group, abx,clearance, dose_level)	

initial_tyc <- tyc_df %>% 
	inner_join(initial_intra_cage, by = c('rows' = 'group')) %>% 
	inner_join(select(initial_intra_cage, group), by = c('columns' = 'group')) %>% 
	mutate(rows = gsub('D.*', '', rows),
		columns = gsub('D.*', '', columns)) %>% 
	separate(rows, c('cage_rows', 'mouse_rows')) %>% 
	separate(columns, c('cage_cols', 'mouse_cols'))
# compare mice against other cages at inital
initial_inter_cage_diversity <- initial_tyc %>% 
	filter(cage_rows != cage_cols) %>% 
	mutate(comparison = 'Initial Between Cages') %>% 
	select(comparison, clearance, abx, dose_level, tyc = distances)
# compare mice within cage at initial
initial_intra_cage_diversity <- initial_tyc %>% 
	filter(cage_rows == cage_cols) %>% 
	mutate(comparison = 'Initial Within Cages') %>% 
	select(comparison, clearance, abx, dose_level, tyc = distances)

compare_time_point <- function(timepoints){
	samples_list <- meta_file %>% 
		filter(time_point %in% timepoints) %>% 
		select(group, clearance, abx, dose_level)	

	output_diversity <- tyc_df %>% 
		inner_join(samples_list, by = c('rows' = 'group')) %>% 
		inner_join(select(samples_list, group), by = c('columns' = 'group')) %>% 
		mutate(rows = gsub('D.*', '', rows),
			columns = gsub('D.*', '', columns)) %>% 
		filter(rows == columns) %>% 
		mutate(comparison = paste(timepoints, collapse = " to ")) %>% 
		select(comparison, clearance, abx, dose_level, tyc = distances)
	return(output_diversity)
}



# compare mice on day 0 to initial
# compare day day 0 to end
# compare intial to end
comparisons <- list(c('Day 0', 'Initial'), c('End', 'Initial'))
tyc_comparisons <- map_dfr(comparisons, ~ compare_time_point(.x))

initial_tyc_df <- rbind(tyc_comparisons, initial_intra_cage_diversity, initial_inter_cage_diversity) 

initial_tyc_df %>% 
	mutate(comparison = factor(comparison, 
		levels = c('Initial Within Cages', 'Initial Between Cages', 
			'Day 0 to Initial', 'Day 0 to End', 'End to Initial'))) %>% 
	filter(clearance != 'Clearing') %>% 
	ggplot(aes(x = comparison, y = tyc, fill = clearance)) +
		geom_boxplot() + 
		theme_bw() + 
		facet_grid(abx~dose_level)

# compare across antibiotic in mice that clear on initial day
# compare across antibiotic in mice that clear on day 0
# compare across antibiotic in mice that clear on end point
cleared_df <- meta_file %>% 
	filter(clearance == 'Cleared', time_point %in% c('End', 'Initial', 'Day 0')) %>% 
	select(group, abx, dose_level, time_point)
across_abx_tyc <- tyc_df %>% 
	inner_join(cleared_df, by = c('rows' = 'group')) %>% 
	inner_join(select(cleared_df, group, abx), by = c('columns' = 'group'),
		suffix = c('_rows','_cols'))
across_abx_tyc %>% 
	filter(abx_rows != abx_cols) %>% 
	filter(time_point == 'Initial')
	ggplot(aes(y = distances, fill = abx_rows)) + 
		geom_boxplot() + 
		theme_bw() + 
		facet_grid(dose_level~time_point)
str(meta_file$time_point)