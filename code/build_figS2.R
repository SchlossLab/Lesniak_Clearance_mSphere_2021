##############
#
# run script to generate Figure S2
#	How does relative abundance differ for the samples with increased alpha diversity?
# 
# Nick Lesniak 02-18-2021
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/sample.final.0.03.subsample.shared
#
##############


library(tidyverse)

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/sample.final.0.03.subsample.shared'

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

# select data from figrue 2B
meta_cef_df <- meta_df %>% 
  filter(abx == 'Cefoperazone',
         clearance != 'Uncolonized' ,
         clearance != 'Clearing',
    time_point %in% c('Initial', 'Day 0', 'End')) %>% 
  select(group, dose, mouse_id, clearance, time_point)

cef_shared <- shared_df %>% 
  filter(Group %in% meta_cef_df$group) %>% 
  pivot_longer(-Group, names_to = 'otu', values_to = 'count') %>% 
  group_by(otu) %>% 
  mutate(present = max(count) > 0) %>% 
  filter(present == T) %>% 
  left_join(meta_cef_df, by = c('Group' = 'group'))

# determine samples with high sobs from fig 2B
increased_div_df <- cef_shared %>% 
  group_by(Group) %>% 
  filter(count > 0) %>% 
  summarise(n_otus = length(count)) %>% 
  mutate(increased_div = n_otus > 120) %>% 
  select(Group, increased_div) %>% 
  mutate(increased_div = factor(increased_div))
                                
levels(increased_div_df$increased_div) <- c(expression(~S[obs]~'< 120'), expression(~S[obs]~'> 120'))

# plot distribution of relative abundance across OTUs
cef_alpha_diff_plot <- cef_shared %>% 
  left_join(increased_div_df) %>% 
  mutate(time_point = factor(time_point, 
                             labels = c('Initial', 'Time~of~Challenge', 'End'),
                             levels = c('Initial', 'Day 0', 'End'))) %>% 
  filter(count > 0) %>% 
  ggplot(aes(x = otu, y = count/15, color = increased_div)) + 
    geom_point(alpha = 0.3) + 
    scale_y_log10() + 
    facet_grid(time_point~increased_div, labeller = label_parsed) + 
    labs(x = 'OTU', y  = 'Relative Abundance') + 
    theme_bw() + 
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

ggsave('results/figures/figure_S2.jpg',
       cef_alpha_diff_plot,
       width = 6, height = 11, units = 'in')
