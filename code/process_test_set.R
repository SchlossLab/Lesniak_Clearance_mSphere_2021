library(tidyr)
library(dplyr)
# From excel file 
# This dataset contains time-series measurements of relative species abundance 
# for the multi-species communities shown in Fig. 3A. The name of each community
# is indicated above each matrix of relative abundances. For example, ER- denotes
# a community containing all species in the synthetic human gut consortium except 
# ER (11 species). "NONE" corresponds to the full 12-member community. The data 
# represents the mean species abundance of each organism 
# (order listed below) at each time point.  {in "Description Sheet"}

test_data <- readxl::read_excel('data/raw/inline-supplementary-material-4-1.xlsx', 
	sheet = 'Dataset EV3', range = 'A1:A104', col_names = F)

# From excel file:
# The relative abundance data is listed in the following order: 
species_labels <- readxl::read_excel('data/raw/inline-supplementary-material-4-1.xlsx', 
	sheet = 'Description', range = 'A4:A15', col_names = F) %>% 
		pull


experiment_label <- test_data[seq(1,nrow(test_data), 8), ] %>% 
	pull
experiment_label <- gsub('\'','', experiment_label)
experiment_label <- gsub('NONE', 'Full', experiment_label)

experiment_data <- test_data[-seq(1,nrow(test_data), 8), ] %>% 
	separate(col = X__1, into = species_labels, sep = '    ') %>% 
	sapply(as.numeric)

clean_df <- data.frame(
	experiment = rep(experiment_label, each = 7),
	time = rep(c(0, 12.4, 25.3, 36, 47.5, 61, 71.4), 13),
	experiment_data,
	stringsAsFactors = F)

write.table(clean_df, 'data/process/mccm_test_data_venturelli_et_al_2018.txt', 
	sep = '\t', quote = FALSE, row.names = FALSE)

## replicate Figure 3A
#clean_df  %>% 
#	gather(species, rel_abun, -experiment, -time) %>% 
#	ggplot(aes(x = time, y = rel_abun, fill = species)) +
#		geom_bar(position = 'stack', stat = 'identity') +
#		facet_wrap(.~experiment, nrow = 2)