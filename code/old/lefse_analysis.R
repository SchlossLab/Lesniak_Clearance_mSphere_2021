library(Hmisc)
library(tidyverse)


# create data frame for lefse 

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
taxonomy_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
taxonomy_function <- 'code/sum_otu_by_taxa.R'
save_dir <- 'results/figures/lefse/'
mothur <- '/mothur/mothur'

# number of samples OTU must be in to be counted
OTU_threshold <- 6

meta_df <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
  mutate(treatment = paste(abx, dose, delayed, cdiff, sep = '_'),
    group = gsub('-', '_', group),
    day = as.numeric(day),
    colonized = ifelse(max_cfu > 10^4, T, F))  %>% 
    filter(abx %in% c('Cefoperazone', 'Streptomycin', 'Clindamycin'),
      clearance != 'Clearing') %>% 
      mutate(time_point = gsub('Day 0', 'Day_0', time_point))
shared_df_raw <- read.table(shared_file, sep = '\t', header = T, stringsAsFactors = F)

source(taxonomy_function)
taxonomy_df <- read.table(taxonomy_file, sep = '\t', header = T, stringsAsFactors = F)
rm(meta_file, shared_file, taxonomy_file, taxonomy_function)

shared_df <- select(shared_df_raw, -label, -numOtus, -Group)
rownames(shared_df) <- shared_df_raw$Group

# remove OTUs present in fewer than specified threshold
samples_w_OTU <- apply(shared_df, 2, function(x) sum(x > 3))
OTUs_above_threshold <- names(samples_w_OTU[samples_w_OTU >= OTU_threshold])
shared_df <- select(shared_df_raw, label, Group, numOtus, 
  one_of(OTUs_above_threshold)) %>% 
  mutate(Group = gsub('-', '_', Group), numOtus = length(OTUs_above_threshold))

# calculate relative abundances
relative_abundance_df <- select(shared_df, -label, -Group, -numOtus)
abundance <- apply(relative_abundance_df, 1, sum)
relative_abundance_df <- data.frame(100 * relative_abundance_df / abundance) %>% 
  mutate(group = shared_df$Group)


rm(samples_w_OTU, OTUs_above_threshold, OTU_threshold, shared_df_raw)

# make sure directory for figures exists
ifelse(!dir.exists(save_dir), 
  dir.create(save_dir), 
  print(paste0(save_dir, ' directory ready')))

# function to run lefse on subset dataframe and produce plot of data
# input_dataframe_name <- 'cef_0.3_clearance_df' # for testing
plot_lefse <- function(input_dataframe_name){
  i <- input_dataframe_name
  current_df <- get(i)
  current_shared <- shared_df %>% 
    filter(Group %in% current_df$group)

  # remove otus that either have 0 or only present in 2 or fewer samples
  present_otus <- current_shared %>% 
    select(-label, -Group, -numOtus) %>% 
    map_dbl(~ sum(. > 0)) %>% 
    which(x = (. > 2)) %>% 
    names
  current_shared <- current_shared %>% 
    select(label, Group, numOtus, one_of(present_otus)) %>% 
    mutate(numOtus = length(present_otus))

  # write files to be used in mothur for lefse analysis
  write_tsv(path = paste0('data/mothur/', i, '.shared'), 
    x = filter(current_shared, Group %in% current_df$group))
  write_tsv(path = paste0('data/mothur/', i, '.design'), 
    x = filter(current_df, group %in% current_shared$Group))

  # run lefse
  system(paste0(mothur, ' "#set.dir(input=data/mothur, output=data/mothur);
    lefse(shared=', i, '.shared, design=', i, '.design);"'))

  # plot lefse results
  lefse_df <- read_tsv(paste0('data/mothur/', i, '.0.03.lefse_summary'))%>% 
    top_n(10, LDA) %>% 
    arrange(desc(LDA)) 

  if(nrow(lefse_df) > 0){
    plot_df <- current_df %>% 
      left_join(relative_abundance_df, by = 'group') %>% 
      select(group, class, one_of(lefse_df$OTU)) %>%  
      gather(OTU, abundance, contains('Otu00')) %>% 
      left_join(select(taxonomy_df, OTU, tax_otu_label), by = 'OTU') %>% 
      left_join(lefse_df, by = 'OTU') %>% 
      mutate(tax_otu_label = gsub('unclassified', 'UC', tax_otu_label),
        tax_otu_label = paste0(tax_otu_label, '\n p = ', pValue))

    lefse_plot <- ggplot(data = plot_df, 
        aes(x = tax_otu_label, y = (abundance + 0.01), color = as.factor(class))) + 
      geom_jitter(position = position_jitterdodge(dodge.width = 0.9), alpha = 0.25) + 
      #geom_boxplot(fill = NA, outlier.color = NA, coef = 0, # whisker length = coef * IQR
      #  position = position_dodge(width = 0.9)) + 
      stat_summary(fun.data = 'median_hilow', geom = 'crossbar', aes(group = class),
        fun.args = (conf.int=0.5), position = position_dodge(width = 0.9)) +
      theme_bw() + labs(x = NULL, y = 'Abundance (counts)', title = paste0('LEfSe on ', i)) + 
      scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
          labels=c("0","0.1","1","10","100")) +
      geom_vline(xintercept=seq(1.5, length(unique(plot_df$tax_otu_label))-0.5, 1), 
               linetype = 'dashed', lwd=0.5, colour="black") + 
      theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) + 
      coord_flip()

    ggsave(paste0(save_dir, 'lefse_plot_', i, '.jpg'), lefse_plot, width = 9, height = 12)
  } else {
    print('No significant LDA')
  }
  return(mutate(lefse_df, analysis = i))
}

# create subsets dataframes with category of interests

# as a whole, what OTUs best explain clearance
cleared_cef_clinda_strep_df <- meta_df %>% 
  filter(cdiff == T, day == 0, colonized == T, 
    clearance %in% c('Cleared', 'Colonized'),
    group %in% shared_df$Group) %>% 
  select(group, class = clearance)

plot_lefse('cleared_cef_clinda_strep_df')

# what are the differences between initial and final in clearing mice for
# clinda - diff between begining and end for mice that clear
clinda_clearance_df <- meta_df %>% 
	filter(abx == 'Clindamycin', clearance == 'Cleared', 
		time_point %in% c('Day_0', 'End')) %>% 
	select(group, class = time_point)
plot_lefse('clinda_clearance_df')

# cef - diff between begining and end for mice that clear
cef_clearance_df <- meta_df %>%
	filter(abx == 'Cefoperazone', clearance == 'Cleared', 
		time_point %in% c('Day_0', 'End')) %>% 
	select(group, class = time_point)
plot_lefse('cef_clearance_df')
# strep - diff between begining and end for mice that clear
strep_clearance_df <- meta_df %>%
	filter(abx == 'Streptomycin', clearance == 'Cleared', 
		time_point %in% c('Day_0', 'End')) %>% 
	select(group, class = time_point)
plot_lefse('strep_clearance_df')


cef_0.3_clearance_initial_df <- meta_df %>% 
  filter(cdiff == T, day == 0, abx == 'Cefoperazone', dose == '0.3', colonized == TRUE) %>% 
  select(group, class = clearance)
# for Cef 0.3, dose with both outcomes, what OTUs best explain clearance
# based on initial community
cef_0.3_clearance_initial_df <- meta_df %>% 
  filter(cdiff == T, day == 0, abx == 'Cefoperazone', dose == '0.3', colonized == TRUE) %>% 
  select(group, class = clearance)
# based on final community
cef_0.3_clearance_end_df <- meta_df %>% 
  filter(cdiff == T, day == last_sample, abx == 'Cefoperazone', dose == '0.3', colonized == TRUE) %>% 
  select(group, class = clearance)
# for all Cef, what OTUs best explain clearance
# based on initial community
cef_all_clearance_initial_df <- meta_df %>% 
  filter(cdiff == T, day == 0, abx == 'Cefoperazone', colonized == TRUE) %>% 
  select(group, class = clearance)
 dose with both outcomes,
# based on final community
cef_all_clearance_end_df <- meta_df %>% 
  filter(cdiff == T, day == last_sample, abx == 'Cefoperazone', colonized == TRUE) %>% 
  select(group, class = clearance)


# for strep 0.5, what OTUs best explain clearance
strep_0.5_clearance_df <- meta_df %>% 
  filter(cdiff == T, day == 0, abx == 'strep', dose == '0.5', colonization == TRUE) %>% 
  select(group, class = clearance)
#  [1] "No significant LDA"
strep_0.5_clearance_at_end_df <- meta_df %>% 
  filter(cdiff == T, day == 9, abx == 'strep', dose == '0.5', colonization == TRUE) %>% 
  select(group, class = clearance)

# for strep 0.1/0.5, what OTUs best explain the difference in antibiotic effect?
strep_0.1_0.5_comparison_df <- meta_df %>% 
  filter(cdiff == T, day == 0, abx == 'strep', dose %in% c('0.1', '0.5')) %>% 
  select(group, class = dose)

# repeat looking at genus level?
# compare community at end point for clearance?
# what explains the variation in abx effect?



# run plotting function using the name of the subsetted dataframe
# expecting two columns, 1 with group - sample names, 1 with class - grouping category
output <- c()
for(i in c('cleared_df', 'colonized_df', 'cef_0.3_colonization_df', 'cef_0.3_clearance_df', 
    'metro_delayed_colonization_df', 'metro_delayed_recovery_df', 'strep_0.5_clearance_df',
    'strep_0.1_0.5_comparison_df', 'strep_0.5_clearance_at_end_df', 'vanc_0.1_clearance_df')){
    output <- rbind(output, plot_lefse(i))
  }
#plot_lefse('cef_0.3_clearance_df')

# create dataframe with the differenced CFU/abundance
delta_df <- meta_df %>% 
  filter(cdiff == T, colonized == T, 
    group %in% shared_df$Group) %>% 
  group_by(mouse_id) %>% 
  arrange(mouse_id , day) %>% 
  select(group, mouse_id, day, CFU, delta_log10cfu) %>% 
  left_join(
    select(shared_df, -numOtus, -label), 
    by = c('group' = 'Group')) %>% 
  mutate_at(vars(contains("Otu0")), function(x) c(NA, diff(x))) %>% 
  filter(day > 0) %>% 
  gather(OTU, delta_abundance, contains('Otu0'))  %>% 
  ungroup
# elimniate otus that do not change
good_otus <- delta_df %>%
  group_by(mouse_id, OTU) %>% 
  summarise(sum_delta_abdnc = sum(abs(delta_abundance))) %>% 
  filter(sum_delta_abdnc != 0) %>% 
  select(mouse_id, OTU) %>% 
  ungroup
# correlate the change in CFU to the change in abundance
delta_corr_otus <- delta_df %>% 
  right_join(good_otus, by = c('mouse_id', 'OTU')) %>% 
  group_by(OTU) %>% 
  nest() %>% 
  # calculate the Spearman correlation value and its P-value for each OTU against
  mutate(comparison = map(data, ~cor.test(.$delta_log10cfu, .$delta_abundance, 
      method = 'spearman')),
    summaries = map(comparison, broom::glance)) %>% 
  unnest(summaries) %>% 
  select(OTU, estimate, p.value) %>% 
  # identify the significant correlations after applying a correction for multiple
  # comparisons
  mutate(p_adjust = p.adjust(p.value, method="BH") < 0.05) %>% 
  filter(p_adjust == TRUE)

abundance_corr_otu <- meta_df %>% 
  filter(cdiff == T, colonized == T, 
    group %in% shared_df$Group) %>% 
  select(group, mouse_id, day, CFU = log10CFU) %>% 
  left_join(
    select(shared_df, -numOtus, -label), 
    by = c('group' = 'Group')) %>% 
  filter(day > 0) %>% 
  gather(OTU, abundance, contains('Otu0'))  %>% 
  group_by(OTU) %>% 
  nest() %>% 
  # calculate the Spearman correlation value and its P-value for each OTU against
  mutate(comparison = map(data, ~cor.test(.$CFU, .$abundance, 
      method = 'spearman')),
    summaries = map(comparison, broom::glance)) %>% 
  unnest(summaries) %>% 
  select(OTU, estimate, p.value) %>% 
  # identify the significant correlations after applying a correction for multiple
  # comparisons
  mutate(p_adjust = p.adjust(p.value, method="BH") < 0.05) %>% 
  filter(p_adjust == TRUE) %>% 
  mutate(est_mag = abs(estimate)) %>% 
  arrange(desc(est_mag))

delta_corr_otus
abundance_corr_otu


#trend_df <- 
meta_df %>% 
  filter(colonized == T, day > 0) %>%
  select(CFU, day, mouse_id) %>% 
  spread(mouse_id, CFU) %>% 
  as.matrix %>% 
  {lm.fit(.[,-1], .[,1])}
    
  #ggplot(aes(x = day, y = log10(CFU + 1))) + geom_point() + coord_cartesian(ylim = c(0,9)) %>% 
    lm.fit(day, CFU) %>% 
    summary(.)['coefficients'] %>% 
    str
    group %in% shared_df$Group) %>% 
  select(group, class = clearance) 



day1 <- c(1,3,1)
day2 <- c(2,2,1)
day3 <- c(3,1,5)
dat <- data.frame(day1,day2,day3)


fits <- lm.fit(cbind(1, seq_len(nrow(dat))), t(dat))
t(coef(fits))

meta_df %>% 
  filter(treatment == 'metro_1_FALSE_TRUE', day == 10)  %>% 
  select(mouse_id, CFU)  




################################################################################
significant_OTUs <- output %>% 
  filter(analysis == 'strep_0.5_clearance_at_end_df') %>% 
  pull(OTU) 

meta_df %>% 
  filter(treatment == 'strep_0.5_FALSE_TRUE') %>% 
  inner_join(select(relative_abundance_df, group, 
      one_of(significant_OTUs))) %>% 
  mutate(Cdiff = log10(CFU)) %>% 
  select(mouse_id, day, Cdiff, one_of(significant_OTUs)) %>% 
  gather(OTU, abundance, Cdiff, one_of(significant_OTUs)) %>% 
  left_join(select(taxonomy_df, OTU, tax_otu_label), by = 'OTU') %>% 
  mutate(tax_otu_label = ifelse(is.na(tax_otu_label), 'C_difficile', tax_otu_label),
    tax_otu_label = gsub('unclassified', 'UC\n', tax_otu_label)) %>% 
  ggplot(aes(x = day, y = abundance, group = mouse_id, color = mouse_id)) + 
    geom_point() + geom_line() +
    facet_grid(tax_otu_label~., scales = 'free_y')

################################################################################
significant_OTUs <- output %>% 
  filter(analysis == 'cef_0.3_clearance_df') %>% 
  pull(OTU) 

meta_df %>% 
  filter(treatment == 'cef_0.3_FALSE_TRUE') %>% 
  inner_join(select(relative_abundance_df, group, 
      one_of(significant_OTUs))) %>% 
  mutate(Cdiff = log10(CFU)) %>% 
  select(mouse_id, day, cage, Cdiff, one_of(significant_OTUs)) %>% 
  gather(OTU, abundance, Cdiff, one_of(significant_OTUs)) %>% 
  left_join(select(taxonomy_df, OTU, tax_otu_label), by = 'OTU') %>% 
  mutate(tax_otu_label = ifelse(is.na(tax_otu_label), 'C_difficile', tax_otu_label),
    tax_otu_label = gsub('unclassified', 'UC', tax_otu_label),
    tax_otu_label = gsub('\\(', '\n\\(', tax_otu_label)) %>% 
  ggplot(aes(x = day, y = abundance, group = mouse_id, color = as.factor(cage))) + 
    geom_point() + geom_line() +
    facet_grid(tax_otu_label~., scales = 'free_y') +
     theme(strip.text.y = element_text(angle = 0))

# Cef 0.3 has a interesting outcome, one set of cages 104 and 105 have a lot of increases in bacteria not obeserved in the other cages with the same treatment and they are highly colonized and do not clear. Where as the other two cages do not have these same increases and lead to either colonization and clearance or no colonization.
# cages 104/105 cages are more similar to those with cef 0.5 treatment
cages_104_105 <- meta_df %>% 
  filter(cage %in% c(104, 105),
    day == 0) %>% 
  select(group, mouse_id) %>% 
  left_join(relative_abundance_df) %>% 
  gather(OTU, abundance, contains('Otu00')) 
OTUs <- cages_104_105 %>% 
  group_by(OTU) %>% 
  summarise(abundance = mean(abundance)) %>% 
  top_n(11, abundance) %>% 
  left_join(select(taxonomy_df, OTU, tax_otu_label), by = 'OTU') %>% 
  select(OTU, tax_otu_label)
# since these two cages seem to be affected differently, theres an increase in Other
# what makes up Other, are there dominant OTUs in Other?
cages_104_105 %>% 
  full_join(OTUs, by = 'OTU') %>% 
  mutate(tax_otu_label = ifelse(is.na(tax_otu_label), 'Other', tax_otu_label)) %>% 
  mutate(tax_otu_label = factor(tax_otu_label, c('Other', OTUs$tax_otu_label[OTUs$tax_otu_label!='Other']))) %>% 
  ggplot(aes(x = mouse_id, y = abundance, fill = tax_otu_label, group = interaction(tax_otu_label,OTU))) +
    geom_bar(stat="identity", position='stack', width = 1, color = "black", size = 0.1) + 
    theme_bw() + labs(x = NULL) + 
    theme(legend.position = 'top', legend.title=element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.text=element_text(size=6)) + 
    scale_fill_brewer(palette = 'Paired')
# doesnt look like any dominant OTUs in Other
#  to eliminate question about dominant in other
#   heat map relative abundances
#   split other bar into individual bars but colored the same

# notes for next steps
#  restore diverse porphyormonadaceae?
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3839982/
#   https://www.frontiersin.org/articles/10.3389/fcimb.2019.00006/full