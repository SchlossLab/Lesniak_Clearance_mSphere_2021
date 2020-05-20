##############
#
# run script to generate plots for Figure 5
#	How do communities predict/classify clearance/colonization?
# 
# Nick Lesniak 04-28-2020
#
#  need files:
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	
#
##############

library(tidyverse)
library(cowplot)
library(pROC)
source("code/R/functions.R")

args <- commandArgs(trailingOnly = TRUE)
model_dir <- as.character(args[1])

# read in data and create key/label to add description to features
meta_file <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_df <- read_tsv(meta_file)[1,] %>% 
	transmute(key = paste(abx, dose_level, 'cage', cage, sep = '_'),
		label = paste0(abx, ' ', dose_level, ' (Cage ', cage, ')'),
		Group = group, mouse = mouse) %>% 
		unique
abx_df <- read_tsv(meta_file) %>% 
	transmute(key = abx, label = abx) %>% 
	unique
tax_file <- 'data/process/abx_cdiff_taxonomy_clean.tsv'
tax_df <- read_tsv(tax_file) %>% 
	mutate(tax_otu_label = gsub('_unclassified', '', tax_otu_label))  %>% 
	select(key = OTU, label = tax_otu_label)
# read in model performance
l2_performance <- paste0('data/process/', model_dir, 
	'/combined_best_hp_results_L2_Logistic_Regression.csv')
l2_performance_df <- read_csv(l2_performance)
samples <- read_csv(paste0('data/process/', model_dir, '_sample_names.txt'))
outcomes <- read_csv(paste0('data/process/', model_dir, '_input_data.csv'))
l2_performance_by_sample <- paste0('data/process/', model_dir, 
	'/combined_all_sample_results_L2_Logistic_Regression.csv')
l2_sample_perf_df <- read_csv(l2_performance_by_sample) %>% 
	inner_join(tibble(Group = samples$Group, cdiff_colonization = outcomes$cdiff_colonization))


# plot the aucs
model_perf_plot <- l2_performance_df %>% 
	rename(`Test\nAUC` = test_aucs, `CV\nAUC` = cv_aucs) %>% 
	pivot_longer(cols = contains('AUC'), values_to = 'AUC', names_to = 'validation', ) %>% 
	ggplot(aes(x = validation, y = AUC, fill = model)) +
		geom_boxplot(alpha=0.5, fatten = 2) +
		geom_hline(yintercept = 0.5, linetype="dashed") +
		coord_flip(ylim = c(0.4, 1)) +
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
		labs(x = NULL, y = "AUROC")


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
    plot <- ggplot(data, aes(reorder(data$label, -data$rank, FUN = median), data$rank)) +
      geom_point(aes(colour= factor(data$sign)), size=1.5) + # datapoints lighter color
      scale_color_manual(values=c("#56B4E9","red3", "#999999"),
      	breaks = c('negative', 'positive', 'zero'),
      	labels = c('Decreased', 'Maintained', 'Null')) +
      stat_summary(fun.y = function(x) median(x), colour = 'black', geom = "point", size = 3) + # Median darker
      coord_flip(ylim=c(0,100)) +
      theme_classic() +
      theme(plot.margin=unit(c(1.5,3,1.5,3),"mm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.title = element_blank(),
          legend.position = c(0.85,0.9),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size = 11, colour='black', face="bold"), 
          axis.text.y=element_text(size = 8, colour='black')) + 
      labs(y = 'Feature ranks')

    return(plot)
}
# -------------------------------------------------------------------->


######################################################################
#--------------Run the functions and plot feature ranks ----------#
######################################################################

label_df <- bind_rows(tax_df, select(meta_df, -mouse, -Group), abx_df)

logit_imp <- get_feature_ranked_files(paste0("data/process/", model_dir, 
	"/combined_L2_Logistic_Regression_feature_ranking.tsv"), 'L2_Logistic_Regression') %>% 
	left_join(label_df, by = c('key'))
logit_graph <- plot_feature_ranks(logit_imp) +
  scale_x_discrete(name = expression(paste(L[2], "-regularized logistic regression"))) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# -------------------------------------------------------------------->


######################################################################
#--------------How does the model do with specific samples ----------#
######################################################################




diff_samples <- c()
for(outcome in c('maintained', 'decreased')){
	# filter samples by true outcome
	data <- l2_sample_perf_df %>% 
		filter(cdiff_colonization %in% outcome)
	# for each mouse, test if prediction for mouse is significantly different
	# from prediction for group
	all_mice <- unique(data$Group)
	for(individual in all_mice){
		subset_sample <- data %>% 
			filter(Group %in% individual) %>% 
			pull(maintained)
		subset_group <- data %>% 
			filter(!Group %in% individual) %>% 
			pull(maintained)
		if(outcome == 'maintained'){
			pvalue <- wilcox.test(subset_sample, subset_group, alternative = 'less')$p.value
		} else {
			pvalue <- wilcox.test(subset_sample, subset_group, alternative = 'greater')$p.value
		}
		diff_samples <- rbind(diff_samples, 
			tibble(sample = individual, pvalue = pvalue))
	}
}

diff_samples <- diff_samples %>% 
	mutate(pvalueBH = p.adjust(pvalue, method = 'BH'),
		significant = case_when(pvalueBH < 0.05 ~ 'Significant',
 			pvalue < 0.05 ~ 'Significant w/o correction',
 			T ~ 'Not Significant'))

medain_df <- l2_sample_perf_df %>% 
	group_by(cdiff_colonization) %>% 
	summarise(maintained = median(maintained)) %>% 
	ungroup

perf_by_sample_plot <- l2_sample_perf_df %>%
	inner_join(diff_samples, by = c('Group' = 'sample')) %>% 
	left_join(meta_df, by = 'Group') %>% 
	mutate(mouse = paste0(' Sample ', mouse, ')'),
		label = str_replace(label, '\\)', mouse)) %>% 
	ggplot(aes(x = Group, y = maintained, color = significant)) +
		facet_wrap(.~cdiff_colonization, scales = 'free') + 
		geom_boxplot() + 
		geom_hline(data = medain_df, aes(yintercept = maintained)) + 
		coord_flip() + 
		scale_color_manual(values = c('grey3', 'red3', 'pink')) + 
		theme_bw() + 
		theme(legend.position = 'none') + 
		labs(x=NULL, y = 'Probabilty Colonization is Maintained',
			caption = 'Red data are significant with multiple comparison correction/nPink data are significant without correction')

# -------------------------------------------------------------------->


######################################################################
#------------- Plot L2 Log Reg Permutation Importance -------------- #
######################################################################

# Read in the cvAUCs, testAUCs for 100 splits as base test_aucs
logit <- read_files(paste0("data/process/", model_dir,
	"/combined_best_hp_results_L2_Logistic_Regression.csv"))
# ----------------------------------------------------------------------------->


# --------  Get the top 20 OTUs that have the largest impact on AUROC ---------->

# Define the function to get the  most important top 20 OTUs
# Order the dataframe from smallest new_auc to largest.
# Because the smallest new_auc means that that OTU decreased AUC a lot when permuted
non_cor_files <-  list.files(path= paste0('data/process/', model_dir), 
  pattern='combined_all_imp_features_non_cor_.*', full.names = TRUE)
for(file_name in non_cor_files){
  importance_data <- read_files(file_name)
  get_interp_info(importance_data, 'L2_Logistic_Regression') %>%
    as.data.frame() %>%
    write_tsv(., paste0("data/process/", model_dir, "/L2_Logistic_Regression_non_cor_importance.tsv"))
}

top_20_otus <-  read_tsv(paste0("data/process/", model_dir, 
	"/L2_Logistic_Regression_non_cor_importance.tsv")) %>%
    arrange(imp) %>%
    head(20)

# Grab the base test auc values for 100 datasplits
data_base <- logit %>%
	select(-cv_aucs) %>%
	mutate(new_auc = test_aucs) %>%
	mutate(names="base_auc") %>%
	select(-test_aucs)

# Have a median base auc value for h-line and for correlated testing
data_base_medians <- logit %>%
	summarise(imp = median(test_aucs), sd_imp = sd(test_aucs)) %>%
	mutate(names="base_auc")

# Get the new test aucs for 100 datasplits for each OTU permuted
# So the full dataset with importance info on non-correlated OTUs
data_full <- read_files(paste0("data/process/", model_dir, 
	"/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv")) %>%
	# Keep OTUs and their AUCs for the ones that are in the top 5 changed
	# (decreased the most)
	filter(names %in% top_20_otus$names) %>%
	inner_join(label_df, by = c('names' = 'key')) %>% 
	mutate(names = label) %>% 
  group_by(names)

# Plot boxplot for the base test_auc values
lowerq <-  quantile(data_base$new_auc)[2]
upperq <-  quantile(data_base$new_auc)[4]
mediandf <-  median(data_base$new_auc) %>%
  data.frame()

# Plot the figure
perm_imp_plot <- ggplot() +
	geom_boxplot(data=data_full, aes(fct_reorder(names, -new_auc), new_auc), alpha=0.7) +
	geom_rect(aes(ymin=lowerq, ymax=upperq, xmin=0, xmax=Inf), fill="gray65") +
	geom_boxplot(data=data_full, aes(x=names, y=new_auc), alpha=0.7) +
	scale_fill_manual(values=cols) +
	geom_hline(yintercept = data_base_medians$imp , linetype="dashed") +
	#geom_hline(yintercept = upperq, alpha=0.5) +
	#geom_hline(yintercept = lowerq, alpha=0.5) +
	coord_flip() +
	theme_classic() +
	labs(y = " AUROC with the OTU permuted randomly") +
	theme(plot.margin=unit(c(1.5,3,1.5,3),"mm"),
		legend.position="none",
		axis.title = element_text(size=10),
		axis.text = element_text(size=10),
		panel.border = element_rect(colour = "black", fill=NA, size=2),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.text.y=element_text(size = 8, colour='black'),
		axis.title.x=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks = element_line(colour = "black", size = 1.1)) + 
	scale_x_discrete(name = "L2_Logistic_Regression") + 
	theme(axis.text.x=element_text(size = 10, colour='black'))




# -------------------------------------------------------------------->


######################################################################
#------------------ Do the features make sense? -------------------- #
######################################################################
# plot abundance by outcome


# -------------------------------------------------------------------->


######################################################################
#-----------------------Save figure -------------------------------- #
######################################################################
#combine with cowplot


ggsave(paste0("results/figures/Figure_5_L2_Logistic_Regression_", model_dir, ".jpg"), 
	plot = plot_grid(
		plot_grid(
			ggdraw(add_sub(model_perf_plot, 'Red point is the post-validation cumulative test AUC', size=7)),
			logit_graph, 
			ggdraw(add_sub(perm_imp_plot, "AUROC with the OTU permuted randomly", size=10, vpadding=grid::unit(0,"lines"), y=5, x=0.65, vjust=4.75)),
			labels = c('A', 'B', 'C'), rel_widths = c(1,2,2), nrow = 1), 
		perf_by_sample_plot, labels = c('', 'C'), ncol = 1),
	width = 14, height = 10, units="in")
