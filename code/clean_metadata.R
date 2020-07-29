###
#
# Clean up issues in metadata file
#	Some decrepencies were discovered in the recording of the data
#
# Dependencies: 
#	* data/raw/abx_cdiff_metadata.tsv
#
# Output: 
#	* data/process/abx_cdiff_metadata_clean.txt
#
###

#load dependencies 
library(dplyr)
library(tidyr)
library(readxl)

#load metadata file
meta_file <- read_xlsx('data/raw/abxD01_IDS.xlsx', sheet = 'all',
		col_types = c('text', 'numeric', 'numeric', 'numeric', 'numeric', 'text', 'text', 'text', 'text', 'text', 'text'),
		na = c('na', 'no')) %>% 
	rename(group = sample, cage = `group#`, mouse = `mouse#`, delayed = time, cdiff = Cdiff, preAbx = preABX) %>% 
	# convert columns to numeric/logical where appropriate
	mutate(dose = as.numeric(gsub('mg/kg', '', dose)),
		delayed = delayed == 'delayed',
		cdiff = cdiff == '630dE',
		preAbx = preAbx == 'preAbx',
		recovDays = as.numeric(gsub('recovD', '', recovDays)))

# some samples have NA and some have 0, 
# many of which are before challenge
# there are 12 samples with NA after day 0 (day 4,5,9, and 10)
na_list <- meta_file %>% 
  filter(is.na(CFU), cdiff == T, day > 0) %>% 
  pull(group)
# "042-3D4"  "048-1D5"  "048-3D5"  "048-4D5"  "050-1D9"  "050-2D9" 
# "050-3D9"  "050-4D9"  "087-4D9"  "098-3D4"  "103-3D10" "111-3D4" 

# visualize where NAs and what the adjacent samples are
# 
#ggsave('~/Desktop/NA_values.jpg',
#  meta_file %>% 
#    filter(group %in% na_list)  %>% 
#    select(cage) %>% unique()  %>% 
#    left_join( select(meta_file, cage, mouse, CFU, day)) %>%
#    mutate(mouse = paste0(cage, '_', mouse), 
#      CFU = CFU + 1,
#      CFU = ifelse(is.na(CFU), 0.1, CFU),
#      shapes = ifelse(CFU < 1, 'NA', 'Actual')) %>% 
#    ggplot(aes(x = day, y = CFU, group = mouse, color = factor(cage))) + 
#      geom_jitter(aes(shape = shapes)) + geom_line() + 
#      scale_y_log10() + 
#      facet_grid(cage~.) +
#      labs(title='Time Series Data from Cages with mice that have an NA of C difficile CFU after Day 0',
#        subtitle = 'Due to log scale, 0s were set to 1 and NAs were set to 0.1. Jitter was added to expose overlayed points')
#)
#
# It seems the NAs are on days that were not plated, but samples were still collected, 
#	as we have 16S data from those samples
# We assume NAs prior to Cdiff challenge had 0 CFUs (Days before 1)
# And we cannot make any assumptions on CFUs post challenge so those samples will remain NA
cfu_cleaned <- meta_file %>% 
	mutate(CFU = case_when(!is.na(CFU) ~ CFU,
		group %in% c(na_list) ~ as.numeric(NA),
		is.na(CFU) ~ 0,
		T ~ as.numeric(NA)))

# few discrepencies with a few mice being mislabeled for preAbx were detected
# so we checked if there are any other discrepencies

# first check if theres any dispencies by cage and day
#discrepency <- cfu_cleaned %>%
#  select(cage, day, abx, dose, delayed, cdiff, preAbx, recovDays) %>%
#  group_by(cage, day) %>% 
#  distinct() 
#discrepency[duplicated(select(discrepency, cage, day)), ]
# three errors
#   group CFU cage mouse day   abx dose delayed cdiff preAbx recovDays
#002-1D-6   0    2     1  -5  none <NA>   FALSE  TRUE   FALSE        NA 
#                                           (should be TRUE preAbx since no abx treatment)
#                                             but this might be noting pbs treatment
#089-1D-6   0   89     1  -6 strep  0.1   FALSE  TRUE   TRUE        NA 
#                       (should be 0.5 dose)
#600-2D-6  NA  600     2  -6   cef  0.5   FALSE  TRUE  FALSE        NA
#                                           (should be TRUE preAbx due to rest of cage)

# check for discrepencies by cage
#sum(cfu_cleaned %>%
#  select(cage, abx, dose, delayed, cdiff) %>%
#  distinct() %>%
#  duplicated())
# no errors detected

# check for any extra days
#sum(cfu_cleaned %>%
#  select(cage, mouse, day) %>%
#  duplicated())
#
# no extra days

# check to make sure group == day/mouse/cage
#check_group <- cfu_cleaned %>%
#  separate(group, c('group','group_day'), sep = 'D') %>%
#  separate(group, c('group_cage','group_mouse'), sep = '-') %>%
#  mutate_at(vars(group_cage,group_mouse, group_day), funs(as.numeric)) %>%
#  select(group_cage, group_mouse, group_day, cage, mouse, day) 
#
#sum(check_group$group_cage != check_group$cage)
#sum(check_group$group_mouse != check_group$mouse)
#sum(check_group$group_day != check_group$day)
# all are equal

# check the variation of cfu by cage to see if theres evidence some were miss entered
#cfu_cleaned %>%
#  select(CFU, cage, day) %>%
#  group_by(cage, day) %>%
#  summarise(sd = sd(CFU)) %>%
#  ggplot(aes(sd)) + 
#    geom_histogram() + 
#    scale_x_log10()
# hard to separate biological variation from random entry errors

cfu_cleaned[cfu_cleaned$group == '600-2D-6', 'preAbx'] <- T
cfu_cleaned[cfu_cleaned$group == '089-1D-6', 'dose'] <- '0.5'  

# create a dataframe with the change in C. difficile CFU by day
cfu_cleaned %<>% 
	mutate(mouse_id = paste(cage, mouse, sep = '_')) %>% 
	group_by(mouse_id) %>% 
	arrange(mouse_id, day) %>% 
	mutate(delta_cfu = CFU - lag(CFU)) %>% 
	ungroup %>% 
	mutate(delta_trend = ifelse(sign(delta_cfu) == -1, -1, 1)) %>% 
	group_by(mouse_id) %>% 
	mutate(last_sample = max(day),
		max_cfu = max(CFU),
		initial_sample = min(day))

# categorize mice based on the colonization levels at the end of the experiment
clearance_df <- cfu_cleaned %>% 
	filter(day == last_sample) %>% 
	select(mouse_id, end_point_cfu = CFU, cdiff, delta_trend, max_cfu) %>% 
	mutate(clearance = case_when(cdiff == F ~ 'Unchallenged',
		max_cfu < 1 ~ 'Uncolonized',
		end_point_cfu == 0 ~ 'Cleared', 
		delta_trend < 0 & end_point_cfu < 100000 ~ 'Clearing', # 10^5 separates the mice 
		T ~ 'Colonized'))
## plot the colonization of the different colonization categories
#meta_file %>% 
#	left_join(clearance) %>% 	
#	ggplot(aes(x = day, y = log10(CFU + 70), color = clearance, group = mouse_id)) + 
#		geom_line(alpha = 0.3) + 
#		facet_wrap(clearance ~ .) + 
#		theme(legend.position="none")
## plot to confirm that we are not missing any mice with negatively trending C. difficile
#meta_file %>% 
#	left_join(clearance) %>% 	
#	filter(clearance == 'colonized') %>% 
#	mutate(Decreasing = ifelse(day %in% c(7:10) & delta_trend < 0, T, F)) %>% 
#	ggplot(aes(x = day, y = log10(CFU + 70), color = Decreasing, group = mouse_id)) + 
#		geom_line(alpha = 0.3) 

cfu_cleaned <- cfu_cleaned %>% 
	left_join(select(clearance_df, mouse_id, clearance), 
		by = 'mouse_id') %>% 
	ungroup %>% 
	mutate(abx = case_when(abx == 'clinda' ~ 'Clindamycin',
			abx == 'strep' ~ 'Streptomycin',
			abx == 'cef' ~ 'Cefoperazone',
			abx == "metro" ~ 'Metronidazole',
			abx == "amp" ~ 'Ampicillin',
			abx == "cipro" ~ 'Ciprofloxacin',
			abx == "none" ~ 'None',
			abx == "vanc" ~ 'Vancomycin'),
		log10CFU = ifelse(CFU == 0, 0, log10(CFU)),
		delta_log10cfu = c(NA, diff(log10CFU)),
		#ifelse(delta_cfu == 0, 0, 
		#	delta_trend * log10(abs(delta_cfu))),
		log10CFU = ifelse(log10CFU == 0, log10(60), log10CFU),
		time_point = factor(case_when(day == 0 ~ 'Day 0',
			day == initial_sample ~ 'Initial',
			day == last_sample ~ 'End',
			T ~ 'Intermediate'), levels = c('Initial', 'Day 0', 'End', 'Intermediate'))) %>% 
	group_by(abx) %>% 
	mutate(dose_level = case_when(dose == min(dose) ~ 'Low',
			dose == max(dose) ~ 'High',
			T ~ 'Mid'),
		group = ifelse(day < 0, paste0(cage, '_', mouse, '_Dminus', abs(day)), paste0(cage, '_', mouse, '_D', day)))

# with cleaned data, filter sample used for this set of experiments
cfu_cleaned <- cfu_cleaned %>% 
	filter(abx %in% c('Clindamycin', 'Cefoperazone', 'Streptomycin')) %>% 
	select(-delayed, -recovDays)

write.table(cfu_cleaned, 'data/process/abx_cdiff_metadata_clean.txt', 
	sep = '\t', quote = FALSE, row.names = FALSE)