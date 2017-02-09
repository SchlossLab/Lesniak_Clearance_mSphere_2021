###################
#
# taxa_labels.R
#
# creates table of tax labels for OTUs
# created by Nicholas Lesniak in 2016
#    
# get taxonomy labels from a vector of OTUs and the taxonomy file
# 
# taxa_level - genus, family, order, class, phylum, kingdom
# otu_subset - vector with OTUs
#    when NULL, defualts to all OTUs
# taxa_file - taxonomy file location
#
###################

get_taxa_labels <- function(taxa_file,  taxa_level='genus', otu_subset=NULL){
  taxa_lvl_num <- which(c('genus','family','order','class','phylum','kingdom')==taxa_level)
  taxa_df <- read.table(taxa_file, header = T, stringsAsFactors=FALSE) # read in taxonomy file
  if(is.null(otu_subset)){ # create list of OTUs
    otu_list <- as.character(taxa_df$OTU) # if no OTUs supplied use all
  } else {
    otu_list <- as.character(otu_subset) # ensure supplied OTUs are read as character
  }
  if (taxa_lvl_num %in% c(1:5)){ # convert taxonomy file list to dataframe
          otu_taxonomy <- taxa_df[taxa_df$OTU %in% otu_list,] # subset taxa_df by OTU list
          otu_taxonomy <- otu_taxonomy[as.numeric(factor(otu_list)), ] # order otu_df to match input
# simplify following with dplyr (suggestion from amanda and marc)
# taxonomy <- mutate(otu_taxonomy,Taxonomy= gsub(otu_taxonomy$Taxonomy, pattern="\(\d*\)",replacement="")) %>% separate(Taxonomy,sep=";", c("phylum", "class", "order", "family", "genus"))
          taxonomy <- sapply(otu_taxonomy$Taxonomy,gsub,pattern="\\(\\d*\\)",replacement="") # remove percentage
          taxonomy <- strsplit(as.character(taxonomy),';',fixed=TRUE) # split taxonomy strings into levels
          taxa_length <- max(c(sapply(taxonomy, length), 6)) # set taxa length to max with a min of 6
          taxonomy <- lapply(taxonomy, `[`, seq_len(taxa_length)) # make all OTUs taxonomy vectors the same length
          taxonomy <- data.frame(do.call('rbind', taxonomy)) # create dataframe rows = otu, col = taxa level
          level <- 7-taxa_lvl_num # set level so input of 1 equals genus column 6
          taxa_out <- as.character(taxonomy[ ,level]) # output selected taxa level
          for (i in level:2){ # replace NA classification with classification from next level 
            next_level <- i-1
            taxa_out[is.na(taxa_out)] <- as.character(taxonomy[is.na(taxa_out),next_level])
          }
          taxa_out <- gsub('_unclassified', '', taxa_out) # remove unclassified from labels
# edit output to have cleaned up OTU and italic formated taxa label
          label_out <- paste0(taxa_out, ' (', 
                              gsub('tu0*', 'TU ', otu_taxonomy$OTU),')') # create labels
          return(data.frame(otu=otu_taxonomy$OTU, taxa= taxa_out, taxa_label = label_out,
                            stringsAsFactors = FALSE))
          } else { # error message if wrong taxa_level input
            print('Error: taxa_class must be genus, family, order, class, phylum')
     }
}
