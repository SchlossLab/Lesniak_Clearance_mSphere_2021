###################
#
# taxa_labels.R
#
# creates table of tax labels for OTUs
# created by Nicholas Lesniak in 2016
#    
# get taxonomy labels from a vector of OTUs and the taxonomy file
# 
# taxa_level - genus, family, order, class, phylum, domain
# otu_subset - vector with OTUs
#    when NULL, defualts to all OTUs
# taxa_file - taxonomy file location
#
###################

library(dplyr)
library(tidyr)

get_taxa_labels <- function(taxa_file,  taxa_level='genus', otu_subset=NULL){
  taxonomic_classification <- c("domain", "phylum", "class", "order", "family", "genus")
  taxa_level_num <- which(taxonomic_classification == taxa_level)
  taxa_df <- read.table(taxa_file, header = T, stringsAsFactors=FALSE) # read in taxonomy file
  if(is.null(otu_subset)){ # create list of OTUs
    otu_list <- as.character(taxa_df$OTU) # if no OTUs supplied use all
  } else {
    otu_list <- as.character(otu_subset) # ensure supplied OTUs are read as character
  }
  if (taxa_level_num %in% c(1:6)){ # convert taxonomy file list to dataframe
    otu_taxonomy <- taxa_df[taxa_df$OTU %in% otu_list,] # subset taxa_df by OTU list
    taxonomy_df <- select(otu_taxonomy, Taxonomy, otu = OTU) %>% 
      separate(Taxonomy,sep="\\(\\d*\\);", taxonomic_classification, extra = 'drop') %>%
      mutate(taxa = .[ , taxa_level_num]) 
    for (i in taxa_level_num:2){ # replace unclassified with next higher level
      next_level <- i-1
      unclassified <- which(taxonomy_df$taxa %in% c('unclassified', 'unclassified_unclassified'))
      taxonomy_df[unclassified, 'taxa'] <- paste0('unclassified_', as.character(taxonomy_df[unclassified, next_level]))
    }
# edit output to have cleaned up OTU and italic formated taxa label
  return(taxonomy_df %>% 
    mutate(tax_otu_label = paste0(taxa, ' (', 
        gsub('tu0*', 'TU ', otu),')'), # create labels
      otu_label = paste0(gsub('tu0*', 'TU ', otu))) # create labels
    )
  } else { # error message if wrong taxa_level input
    print('Error: taxa_class must be genus, family, order, class, phylum, or domain')
  }
}
