###################
#
# sum_otu_by_taxa.R
#
# created by Nicholas Lesniak in 2016
#    
# function outputs a dataframe with OTUs summed by taxa group per sample
# 
# taxonomy_file - taxonomy file location
# otu_df - dataframe with samples as rownames and otus in columns
# taxa_level - genus, family, order, class, phylum, kingdom
#
###################
source('../code/taxa_labels.R') # source function to get taxonomy of OTUs

sum_otu_by_taxa <- function(taxonomy_file, otu_df = shared_file, taxa_level = 'genus'){
  taxa_df <- get_taxa_labels(taxa_file = taxonomy_file, taxa_level = taxa_level, 
                             otu_subset = colnames(otu_df)[grep('Otu\\d', colnames(otu_df))])

  taxa_groups <- as.character(unique(taxa_df$taxa)) # get taxonomy groups
# simplify following with dplyr (following is pat/matts suggestion)
  output_dataframe <- otu_df %>% 
    mutate(group = rownames(otu_df)) %>% 
    gather(otu, counts, contains('Otu0')) %>% 
    inner_join(select(taxa_df, otu, taxa), by = 'otu') %>% 
    group_by(group, taxa) %>%
    summarize(total=sum(counts)) %>% 
    spread(taxa, total) %>% 
    ungroup
  return(output_dataframe)
  }