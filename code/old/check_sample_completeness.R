# some mice are missing a few samples but when combining with shared file, 
# many samples are being lost. let's see where the samples are being lost
library(tidyverse)
# read in the meta and shared files
meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
	select(group, cage, mouse, day, CFU, cdiff, abx, dose, delayed) %>% 
	unite(treatment, abx, dose, delayed) %>% 
	filter(cdiff == T, day >= 0, treatment != 'none_NA_FALSE')
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
shared_file <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F)
source('code/sum_otu_by_taxa.R')
taxonomy_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy'
shared_by_genus <- sum_otu_by_taxa(taxonomy_file = taxonomy_file, 
	otu_df = shared_file, 
	taxa_level = 'genus')
meta_shared <- meta_file %>% 
	mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
	inner_join(shared_by_genus, by = c('group' = "Group")) %>% 
	select(-group)
meta_shared <- meta_shared %>% 
	full_join(data.frame(unique(select(meta_shared, treatment, unique_id)), day = 0), by = c('unique_id', 'day', 'treatment'))


# how many days are missing from each treatment (from the meta file)
days_df <- c()
for(i in unique(meta_file$treatment)){
	total_mice <- length(meta_file %>% 
		mutate(unique_id = paste(mouse, cage, sep = '_')) %>% 
		filter(treatment == i) %>% 
		pull(unique_id) %>% 
		unique)
	day_counts <- meta_file %>%
			filter(treatment == i) %>%
			count(day) %>%
			pull(n) %>% t
	temp <- data.frame(i, total_mice, total_mice - day_counts)
	colnames(temp) <- c('treatment', 'total mice', as.character(0:10))
	days_df <- rbind(days_df, temp)
}
days_df$meta_samples_missing <- apply(days_df[ , 4:13], 1, sum)
#          treatment total mice 0 1 2 3 4 5 6 7 8 9 10 meta_samples_missing
#    clinda_10_FALSE         11 0 0 0 0 1 0 0 0 0 0  0                    1
#   vanc_0.625_FALSE          6 0 0 0 1 1 0 0 0 0 0  0                    2
#     cipro_10_FALSE          5 0 0 0 0 0 1 0 0 0 0  0                    1
#      amp_0.5_FALSE          9 2 0 0 0 1 0 0 0 0 4  0                    5
#      cef_0.5_FALSE          6 0 0 0 0 0 0 0 0 0 0  0                    0
#      metro_1_FALSE          7 0 0 3 0 0 0 0 0 0 1  1                    5
#      cef_0.3_FALSE         13 0 0 0 0 1 0 0 1 0 0  0                    2
#      cef_0.1_FALSE          6 0 0 0 0 0 0 0 0 0 0  0                    0
#    strep_0.1_FALSE         10 0 0 0 0 0 0 0 0 0 0  0                    0
#      strep_5_FALSE          8 0 0 1 0 0 0 0 0 0 0  0                    1
#    strep_0.5_FALSE          9 0 0 0 0 0 4 0 1 0 0  0                    5
#     vanc_0.3_FALSE          8 0 0 1 3 0 0 0 0 0 0  0                    4
#     vanc_0.1_FALSE          9 0 0 2 1 0 0 0 0 0 1  0                    4
#       metro_1_TRUE         14 0 0 0 0 0 0 0 0 0 0  0                    0
#       amp_0.5_TRUE         14 1 0 0 0 1 0 0 0 0 0  0                    1

# how many samples go missing when joining with the shared file
days_df_sh <- c()
for(i in unique(meta_shared$treatment)){
	total_mice <- length(meta_shared %>% 
		mutate(unique_id = paste(mouse, cage, sep = '_')) %>% 
		filter(treatment == i) %>% 
		pull(unique_id) %>% 
		unique)
	day_counts <- meta_shared %>%
			filter(treatment == i) %>%
			count(day) %>%
			pull(n) %>% t
	temp <- data.frame(i, total_mice, total_mice - day_counts)
	colnames(temp) <- c('treatment', 'total mice', as.character(0:10))
	days_df_sh <- rbind(days_df_sh, temp)
}
days_df_sh$shared_samples_missing <- apply(days_df_sh[ , 4:13], 1, sum)
#          treatment total mice 0 1 2 3 4 5 6 7 8 9 10 shared_samples_missing
#    clinda_10_FALSE         11 0 0 0 0 1 0 0 0 0 0  0                      1
#   vanc_0.625_FALSE          6 0 0 0 1 2 0 0 0 0 0  0                      3
#     cipro_10_FALSE          5 0 0 0 0 0 1 0 0 0 0  0                      1
#      amp_0.5_FALSE          9 4 2 0 2 3 2 3 2 0 5  0                     19
#      cef_0.5_FALSE          6 4 2 3 1 0 1 1 1 2 3  0                     14
#      metro_1_FALSE          7 0 1 4 2 0 0 1 0 0 1  1                     10
#      cef_0.3_FALSE         13 1 0 0 0 1 0 0 1 0 0  0                      2
#      cef_0.1_FALSE          6 0 0 0 0 0 0 0 0 0 0  0                      0
#    strep_0.1_FALSE         10 0 1 3 3 2 2 1 0 2 0  1                     15
#      strep_5_FALSE          8 0 0 2 0 0 0 1 0 0 0  0                      3
#    strep_0.5_FALSE          9 0 1 1 0 1 6 0 3 1 1  1                     15
#     vanc_0.3_FALSE          8 0 0 1 3 0 0 1 0 4 2  1                     12
#     vanc_0.1_FALSE          9 0 0 2 2 0 0 1 1 1 1  1                      9
#       metro_1_TRUE         14 1 0 1 2 3 0 0 2 1 2  2                     13
#       amp_0.5_TRUE         14 1 6 0 1 5 3 1 3 0 0  1                     20
days_df_sh$lost <- days_df_sh$shared_samples_missing - days_df$meta_samples_missing
# samples lost to sequencing are
#          treatment lost
#    clinda_10_FALSE    0
#   vanc_0.625_FALSE    1
#     cipro_10_FALSE    0
#      amp_0.5_FALSE   14
#      cef_0.5_FALSE   14
#      metro_1_FALSE    5
#      cef_0.3_FALSE    0
#      cef_0.1_FALSE    0
#    strep_0.1_FALSE   15
#      strep_5_FALSE    2
#    strep_0.5_FALSE   10
#     vanc_0.3_FALSE    8
#     vanc_0.1_FALSE    5
#       metro_1_TRUE   13
#       amp_0.5_TRUE   19

# what is the distribution of lost days?
apply(days_df[,3:13], 2, sum)
#day 		0  1  2  3  4  5  6  7  8  9 10 
#missing 	3  0  7  5  5  5  0  2  0  6  1 
apply(days_df_sh[,3:13], 2, sum)
#day 		0  1  2  3  4  5  6  7  8  9  10 
#missing	11 13 17 17 18 15 10 13 11 15  8 
# lost to sequencing by day
#day 		0  1  2  3  4  5  6  7  8  9 10 
#missing 	8 13 10 12 13 10 10 11 11  9  7 
# day 0 is not included since it is used as the separator

meta_shared	%>%
	group_by(unique_id) %>%
	count(unique_id) %>%
	filter(n == 11) %>%
	nrow

# theres 135 mice in total but only 61 mice have data for all days
meta_shared	%>%
	group_by(unique_id) %>%
	count(unique_id) %>% ungroup() %>%
	full_join(meta_shared) %>%
	filter(n == 11) %>% 
	select(unique_id, treatment) %>%
	unique %>%
	count(treatment) %>%
	select(treatment, n_mice_w_all_days = n) %>%
	full_join(select(days_df, treatment, `total mice`, meta_samples_missing)) %>%
	full_join(select(days_df_sh, treatment, lost))
# each treatment has the following number of mice with all days
# treatment        n_mice_w_all_days  total mice  meta_samples_missing  lost
# amp_0.5_TRUE                     3           14                    1    19
# cef_0.1_FALSE                    6            6                    0     0
# cef_0.3_FALSE                   11           13                    2     0
# cef_0.5_FALSE                    2            6                    0    14
# cipro_10_FALSE                   4            5                    1     0
# clinda_10_FALSE                 10           11                    1     0
# metro_1_FALSE                    3            7                    5     5
# metro_1_TRUE                     5           14                    0    13
# strep_0.1_FALSE                  6           10                    0    15
# strep_5_FALSE                    5            8                    1     2
# vanc_0.1_FALSE                   4            9                    4     5
# vanc_0.3_FALSE                   2            8                    4     8
# vanc_0.625_FALSE                 3            6                    2     1
# amp_0.5_FALSE                   NA            9                    5    14
# strep_0.5_FALSE                 NA            9                    5    10

days_by_mouse <- c()
for(i in unique(meta_shared$treatment)){
	total_mice <- length(meta_shared %>% 
		mutate(unique_id = paste(mouse, cage, sep = '_')) %>% 
		filter(treatment == i) %>% 
		pull(unique_id) %>% 
		unique)
	day_counts <- meta_shared %>%
			filter(treatment == i) %>%
			count(day) %>%
			pull(n) %>% t
	temp <- data.frame(i, total_mice, total_mice - day_counts)
	colnames(temp) <- c('treatment', 'total mice', as.character(0:10))
	days_df_sh <- rbind(days_df_sh, temp)
}


mouse_by_day <- select(meta_shared, treatment, unique_id)  %>% 
	unique %>% 
	count(treatment) %>% 
	rename(n_mice = n)
for(i in 0:10){
	new_day <- meta_shared %>% 
		filter(day %in% 0:i) %>% 
		group_by(treatment) %>% 
		count(unique_id) %>% 
		filter(n == (i + 1)) %>% 
		count(treatment)
	colnames(new_day)[2] <- paste0('day_', i)
	mouse_by_day <- full_join(mouse_by_day, new_day, by = 'treatment')
}
mouse_by_day[is.na(mouse_by_day)] <- 0
mouse_by_day <- rbind(mouse_by_day, c('total', apply(mouse_by_day[,-1], 2, sum)))	
ggsave('scratch/mouse_by_day.jpg',
	mouse_by_day %>% 
		filter(treatment != 'total') %>% 
		gather(day, mice, -treatment, -n_mice) %>% 
		separate(day, c('x', 'day')) %>%
		mutate(mice = as.numeric(mice), day = as.numeric(day)) %>%  
		#mutate(mice = as.numeric(mice)/as.numeric(n_mice), day = as.numeric(day)) %>%  
		ggplot(aes(x = day, y = mice, color = treatment)) + 
			geom_line() + 
			theme_bw() + 
			labs(title = 'Number of mice per treatment with all previous days')
			)

treatment			n_mice	day_0	day_1	day_2	day_3	day_4	day_5	day_6	day_7	day_8	day_9	day_10
amp_0.5_FALSE		9		9		7		7		6		4		4		3		3		3		0		0     
amp_0.5_TRUE		14		14		8		8		7		6		4		3		3		3		3		3     
cef_0.1_FALSE		6		6		6		6		6		6		6		6		6		6		6		6     
cef_0.3_FALSE		13		13		13		13		13		12		12		12		11		11		11		11    
cef_0.5_FALSE		6		6		4		3		3		3		3		2		2		2		2		2     
cipro_10_FALSE		5		5		5		5		5		5		4		4		4		4		4		4     
clinda_10_FALSE		11		11		11		11		11		10		10		10		10		10		10		10    
metro_1_FALSE		7		7		6		3		3		3		3		3		3		3		3		3     
metro_1_TRUE		14		14		14		13		11		8		8		8		6		5		5		5     
strep_0.1_FALSE		10		10		9		7		7		7		6		6		6		6		6		6     
strep_0.5_FALSE		9		9		8		7		7		6		0		0		0		0		0		0     
strep_5_FALSE		8		8		8		6		6		6		6		5		5		5		5		5     
vanc_0.1_FALSE		9		9		9		7		5		5		5		5		4		4		4		4     
vanc_0.3_FALSE		8		8		8		7		4		4		4		3		3		2		2		2     
vanc_0.625_FALSE	6		6		6		6		5		3		3		3		3		3		3		3     
total				135		135		122		109		99		88		78		73		69		67		64		64    

#treatment_subset <- 'cef_0.1_FALSE'
#treatment_subset <- 'clinda_10_FALSE'