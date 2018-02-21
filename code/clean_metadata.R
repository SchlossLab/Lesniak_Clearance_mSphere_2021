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

#load metadata file
meta_file   <- 'data/raw/abx_cdiff_metadata.tsv'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F)


# convert all dosages to same format
meta_file$dose[meta_file$dose == '10mg/kg'] <- 10
meta_file$dose <- as.numeric(meta_file$dose)

# change vanc so all are same format
meta_file$abx[meta_file$abx == 'vanc '] <- 'vanc'


# some samples have NA and some have 0, 
# many of which are before challenge
# there are 12 samples with NA after day 0 (day 4,5,9, and 10)
na_list <- c(meta_file %>% 
  filter(is.na(CFU), cdiff == T, day > 0) %>% 
  pull(group))
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
# And we cannot make any assumptions on CFUs post challenge so those samples will be removed
cfu_cleaned <- meta_file %>% 
  filter(!group %in% c(na_list)) %>% 
  replace_na(list(CFU = 0))

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


write.table(cfu_cleaned, 'data/process/abx_cdiff_metadata_clean.txt', 
	sep = '\t', quote = FALSE, row.names = FALSE)