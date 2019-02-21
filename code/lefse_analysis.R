
library(tidyverse)
#library(cowplot)
#library(RColorBrewer)
#library(bsselectR)


# create data frame for lefse 

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
taxonomy_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
taxonomy_function <- 'code/sum_otu_by_taxa.R'

# number of samples OTU must be in to be counted
OTU_threshold <- 6

meta_df <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F) %>% 
  mutate(treatment = paste(abx, dose, delayed, cdiff, sep = '_'),
    mouse_id = paste(cage, mouse, sep = '_'),
    group = gsub('-', '_', group)) 
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
rm(samples_w_OTU, OTUs_above_threshold, OTU_threshold, shared_df_raw)

# create category of colonize (max CFU > 10^4) and cleared (min of colonized < 10^4)
colonized_mice <- meta_df %>% 
  group_by(mouse_id) %>% 
  summarise(max_cfu = max(CFU)) %>% 
  mutate(colonization = ifelse(max_cfu > 10^4, T, F))

colonized_df <- meta_df %>% 
  left_join(colonized_mice, by = 'mouse_id') %>% 
  filter(cdiff == T, day == 0, 
    group %in% shared_df$Group) %>%
  select(group, colonization)

write_tsv(path = 'data/mothur/colonized.shared', 
  x = filter(shared_df, Group %in% colonized_df$group))
write_tsv(path = 'data/mothur/colonized.design', x = colonized_df)

cleared_mice <- meta_df %>% 
  right_join(
    filter(colonized_mice, colonization == T),
    by = 'mouse_id') %>% 
  filter(day > 0) %>% 
  group_by(mouse_id) %>% 
  summarise(min_cfu = min(CFU)) %>% 
  mutate(clearance = ifelse(min_cfu < 10^3, T, F))

cleared_df <- meta_df %>% 
  left_join(colonized_mice, by = 'mouse_id') %>% 
  left_join(cleared_mice, by = 'mouse_id') %>% 
  filter(cdiff == T, colonization == T, day == 0, 
    group %in% shared_df$Group) %>% 
  select(group, clearance) 

write_tsv(path = 'data/mothur/cleared.shared', 
  x = filter(shared_df, Group %in% cleared_df$group))
write_tsv(path = 'data/mothur/cleared.design', x = cleared_df)


# run lefse

system('/mothur/mothur "#set.dir(input=data/mothur, output=data/mothur);
  lefse(shared=colonized.shared, design=colonized.design);
  metastats(shared=colonized.shared, design=colonized.design);
  lefse(shared=cleared.shared, design=cleared.design);
  metastats(shared=cleared.shared, design=cleared.design);"')

# plot lefse results

lefse_df <- read_tsv('data/mothur/colonized.0.03.lefse_summary')
lefse_df <- lefse_df  %>% 
  top_n(20, LDA) %>% 
  arrange(desc(LDA)) 

library(ggforce)

colonized_df %>% 
  left_join(shared_df, by = c('group' = 'Group')) %>% 
  select(group, colonization, one_of(lefse_df$OTU)) %>%  
  gather(OTU, abundance, contains('Otu00')) %>% 
  left_join(select(taxonomy_df, OTU, tax_otu_label), by = 'OTU') %>% 
  left_join(lefse_df, by = 'OTU') %>% 
  ggplot(aes(x = colonization, y = abundance, color = colonization)) + 
    geom_sina() +
    stat_summary(fun.y = 'median', color = 'black', geom = 'point') +
    facet_wrap(tax_otu_label~., scales = 'free_y') + 
    theme_bw()
  head
  




