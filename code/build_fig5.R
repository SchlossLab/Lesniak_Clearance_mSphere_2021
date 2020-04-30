##############
#
# run script to generate plots for Figure 5
#	How do communities predict/classify clearance/colonization?
# 
# Nick Lesniak 04-28-2020
#
#  need files:
#	data/process/abx_cdiff_metadata_clean.txt
#	data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	code/sum_otu_by_taxa.R
#
##############

library(tidyverse)
library(cowplot)
source("code/learning/functions.R")


meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared'
subsampled_shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
l2_performance <- 'data/process/combined_best_hp_results_L2_Logistic_Regression.csv'


abx_color <- tibble(abx = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
	color = c('#D37A1F', '#3A9CBC', '#A40019'))
# read in data
tax_df <- read_tsv(tax_file)
network_labels <- tax_df %>% 
	mutate(otu_number = gsub('OTU ', '', otu_label),
		tax_otu_label = gsub('_unclassified', '', tax_otu_label),
		tax_otu_label = gsub(' \\(', '\\\n\\(', tax_otu_label)) %>% 
	pull(tax_otu_label)
meta_df   <- read_tsv(meta_file) %>% 
	filter(abx %in% c('Clindamycin', 'Streptomycin', 'Cefoperazone'),
		cdiff == T, clearance != 'uncolonized', day > 0,)
shared_df <- read_tsv(subsampled_shared_file) %>% 
	select(-numOtus, -label) %>% 
	filter(Group %in% meta_df$group)
l2_performance_df <- read_csv(l2_performance)


model_perf_plot <- l2_performance_df %>% 
	pivot_longer(cols = contains('_aucs'), values_to = 'AUC', names_to = 'validation', ) %>% 
	ggplot(aes(x = validation, y = AUC, fill = model)) +
		geom_boxplot(alpha=0.5, fatten = 2) +
		geom_hline(yintercept = 0.5, linetype="dashed") +
		coord_flip() +
		scale_y_continuous(name = "AUROC",
		                   breaks = seq(0.4, 1, 0.1),
		                   limits=c(0.4, 1),
		                   expand=c(0,0)) +
		theme_bw() +
		theme(plot.margin=unit(c(0,1.1,0,0),"cm"),
		      legend.justification=c(1,0),
		      legend.position='none',
		      panel.grid.major.y = element_blank(),
		      panel.grid.major.x = element_line( size=0.6),
		      panel.grid.minor = element_blank(),
		      panel.background = element_blank(),
		      text = element_text(size = 12),
		      axis.text.x=element_text(size = 10, colour='black'),
		      axis.text.y=element_text(size = 10, colour='black'),
		      axis.title.y=element_text(size = 10),
		      axis.title.x=element_text(size = 10),
		      panel.border = element_rect(linetype="solid", colour = "black", fill=NA, size=1.5)) + 
		labs(x = NULL)


## Feature importance

# ------------------- Re-organize feature importance  ----------------->
# This function:
#     1. Takes in a combined (100 split) feature rankings for each model) and the model name
#     2. Returns the top 5 ranked (1-5 with 1 being the highest ranked) OTUs (ranks of the OTU for 100 splits)
get_feature_ranked_files <- function(file_name, model_name){
  importance_data <- read_tsv(file_name)
  ranks <- get_interp_info(importance_data, model_name) %>%
    as.data.frame()
  return(ranks)
}

# This function:
#     1. Top 5 ranked (1-5 lowest rank) OTUs (ranks of the OTU for 100 splits)
#     2. Returns a plot. Each datapoint is the rank of the OTU at one datasplit.
                        
plot_feature_ranks <- function(data){
    # Plot from highest median ranked OTU to least (only top 5) and thir ranks that lay between 1-100
    # Rank 1 is the highest rank
    plot <- ggplot(data, aes(reorder(data$key, -data$rank, FUN = median), data$rank)) +
      geom_point(aes(colour= factor(data$sign)), size=1.5) + # datapoints lighter color
      scale_color_manual(values=c("#56B4E9","red3", "#999999")) +
      stat_summary(fun.y = function(x) median(x), colour = 'black', geom = "point", size = 3) + # Median darker
      coord_flip(ylim=c(0,100)) +
      theme_classic() +
      theme(plot.margin=unit(c(1.5,3,1.5,3),"mm"),
          legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size = 11, colour='black', face="bold"), 
          axis.text.y=element_text(size = 8, colour='black'))
    return(plot)
}

######################################################################
#--------------Run the functions and plot feature ranks ----------#
######################################################################


logit_imp <- get_feature_ranked_files("data/process/combined_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression") %>% 
	inner_join(select(tax_df, OTU, tax_otu_label), by = c('key' = 'OTU')) %>% 
	select(-key) %>% 
	rename(key = tax_otu_label)
logit_graph <- plot_feature_ranks(logit_imp) +
  scale_x_discrete(name = expression(paste(L[2], "-regularized logistic regression"))) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# -------------------------------------------------------------------->


######################################################################
#-----------------------Save figure as .pdf ------------------------ #
######################################################################
#combine with cowplot

linear <- plot_grid(, labels = c("A"), align = 'v', ncol = 1)

ggsave("results/figures/Figure_5.jpg", 
	plot = plot_grid(
		model_perf_plot,
		ggdraw(add_sub(logit_graph, "Feature ranks", vpadding=grid::unit(0,"lines"), y=5, x=0.7, vjust=4.75, size=15)),
		labels = c('A', 'B'), rel_widths = c(1,2)), 
	width = 9, height = 5, units="in")
