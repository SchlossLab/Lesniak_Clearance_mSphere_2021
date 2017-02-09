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

sum_otu_by_taxa <- function(taxonomy_file, otu_df, taxa_level = 'genus'){
  taxa_df <- get_taxa_labels(taxa_file = taxonomy_file, taxa_level = taxa_level, 
                             otu_subset = colnames(otu_df)[grep('Otu\\d', colnames(otu_df))])
# change to read in with stringasfactor = F and eliminate as.character
  taxa_groups <- as.character(unique(taxa_df$taxa)) # get taxonomy groups
# simplify following with dplyr (following is pat/matts suggestion)
#  otu_df$sample <- rownames(otu_df)
#  output_dataframe <- otu_df %>%
#    select(-label, -numOtus) %>%
#    gather("otu", "counts", 1:nrow(taxa_df)) %>%
#    inner_join(.,taxa_df, by=c("otu"="otu")) %>%
#    group_by(sample, taxa) %>%
#    summarize(total=sum(counts))
  output_dataframe <- sapply(taxa_groups, function(i){
    group <- as.character(taxa_df$otu[taxa_df$taxa %in% i]) # get otus of each group
    if(length(group) < 2){ 
      otu_df[,group] # if only one otu in group, output column
      } else {
        apply(otu_df[,group],1,sum) # sum subset columns
        }
    })
  output_dataframe <- cbind(sample = rownames(output_dataframe),
                            data.frame(output_dataframe))
  return(output_dataframe)
  }