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
library(ggtext)
source("code/R/functions.R")

#args <- commandArgs(trailingOnly = TRUE)
model_dir <- 'otu'

# read in data and create key/label to add description to features
meta_file <- 'data/process/abx_cdiff_metadata_clean.txt'
meta_df <- read_tsv(meta_file) %>% 
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
# read in model performance - from model 50/50 split not by cage
l2_performance <- paste0('data/process/', model_dir, 
	'/combined_best_hp_results_L2_Logistic_Regression.csv')
l2_performance_df <- read_csv(l2_performance)
#samples <- read_csv(paste0('data/process/', model_dir, '_sample_names.txt'))
#outcomes <- read_csv(paste0('data/process/', model_dir, '_input_data.csv'))
l2_performance_by_sample <- paste0('data/process/', model_dir, 
	'/combined_all_sample_results_L2_Logistic_Regression.csv')
l2_sample_perf_df <- read_csv(l2_performance_by_sample) #%>% 
	#inner_join(tibble(Group = samples$Group, cdiff_colonization = outcomes$cdiff_colonization))


# plot the aucs
model_perf_plot <- l2_performance_df %>% 
	rename(`Test\nAUC` = test_aucs, `CV\nAUC` = cv_aucs) %>% 
	pivot_longer(cols = contains('AUC'), values_to = 'AUC', names_to = 'validation', ) %>% 
	ggplot(aes(x = validation, y = AUC)) +
		geom_boxplot(alpha=0.5, fatten = 2) +
		geom_hline(yintercept = 0.5, linetype="dashed") +
		coord_cartesian(ylim = c(0.4, 1)) +
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

label_df <- bind_rows(tax_df, select(meta_df, -mouse, -Group), abx_df) %>% 
	mutate(label = paste0('*', label),
		label = gsub(' \\(', '* \\(', label))

logit_imp <- get_feature_ranked_files(paste0("data/process/", model_dir, 
	"/combined_L2_Logistic_Regression_feature_ranking.tsv"), 'L2_Logistic_Regression') %>% 
	left_join(label_df, by = c('key'))
logit_graph <- plot_feature_ranks(logit_imp) +
  scale_x_discrete(name = expression(paste(L[2], "-regularized logistic regression"))) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
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
top_features_file <- paste0("data/process/", model_dir, 
	"/L2_Logistic_Regression_non_cor_importance.tsv")
if(!file.exists(top_features_file)){
	importance_data <- read_files(paste0("data/process/", model_dir,
		"/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv"))
	get_interp_info(importance_data, 'L2_Logistic_Regression') %>%
	    as.data.frame() %>%
	    write_tsv(., paste0("data/process/", model_dir, "/L2_Logistic_Regression_non_cor_importance.tsv"))
	}

top_features <-  read_tsv(top_features_file) %>%
    arrange(imp) %>% 
    head(36)

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
data_full <- read_csv(paste0("data/process/", model_dir, 
	"/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv")) %>%
	# Keep OTUs and their AUCs for the ones that are in the top 5 changed
	# (decreased the most)
	filter(names %in% top_features$names) %>%
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
	labs(y = "AUROC with the OTU permuted randomly",
		x = expression(paste(L[2], "-regularized logistic regression"))) +
	theme(plot.margin=unit(c(1.5,3,1.5,3),"mm"),
		legend.position="none",
		axis.title = element_text(size=10),
		axis.text = element_text(size=10),
		panel.border = element_rect(colour = "black", fill=NA, size=2),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.text.y=element_markdown(),
		axis.ticks = element_line(colour = "black", size = 1.1)) +
	theme(axis.text.x=element_text(size = 10, colour='black'))

# -------------------------------------------------------------------->


######################################################################
#------------------ Do the features make sense? -------------------- #
######################################################################
# plot abundance by outcome
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
full_meta_df <- read_tsv(meta_file)
shared_df <- read_tsv(shared_file) 
total_abundance <- sum(shared_df[1,-c(1:3)])
# create dataframe with top otus, true outcome, abx
top_otu_abundance <- shared_df %>% 
	inner_join(distinct(select(l2_sample_perf_df, Group, clearance)), by = 'Group') %>% 
	inner_join(select(full_meta_df, Group = group, abx), by = 'Group') %>% 
	select(Group, clearance, abx, one_of(top_features$names)) %>% 
	pivot_longer(cols = c(-Group, -clearance, -abx), names_to = 'otu', values_to = 'abundance') %>% 
	mutate(abundance = (100 * abundance/total_abundance) + .05) %>% 
	inner_join(label_df, by = c('otu' = 'key')) # add column with otu label

# plot otus by decreasing importance and distribution of abundance by abx/outcome
top_otu_abundance_plot <- top_otu_abundance %>% 
	inner_join(data_full, by = 'label') %>% # add perm importance aucs to order otus
	ggplot(aes(x = fct_reorder(names, -new_auc), y = abundance, color = clearance)) + 
		geom_boxplot() + 
		scale_y_log10(limits = c(0.04,100),
		   		breaks = c(0.01, 0.1, 1, 10, 100),
		   		labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) + 
		coord_flip() + 
		facet_wrap(.~abx) + 
		theme_bw() +
		labs(y = 'Percent Relative Abundance', x = NULL, color = NULL) + 
		theme(legend.position = 'bottom',
			panel.grid.minor.y = element_blank(),
			axis.text.y = element_markdown(),
			axis.text.x = element_markdown())

# -------------------------------------------------------------------->

######################################################################
#--------------How does the model do with specific samples ----------#
######################################################################
# plot the distribution of probability of clearance for each sample
perf_by_sample_plot <- l2_sample_perf_df %>%
	left_join(meta_df, by = 'Group') %>% 
	mutate(mouse = paste0(' Sample ', mouse, ')'),
		label = str_replace(label, '\\)', mouse)) %>% 
	group_by(label) %>% 
	mutate(median_p = median(Cleared)) %>% 
	ggplot(aes(x = reorder(label, -median_p), y = Cleared, color = clearance)) +
		geom_boxplot() + 
		coord_flip() + 
		theme_bw() + 
		theme(legend.position = c(0.2, 0.075),
			panel.grid.minor.x = element_blank()) + 
		labs(x=NULL, y = 'Probabilty Colonization is Cleared', color = NULL)
# select samples that are inbetween outcomes, ones most likely to be misclassifies
misclass_samples <- l2_sample_perf_df %>%
	left_join(meta_df, by = 'Group') %>% 
	mutate(mouse = paste0(' Sample ', mouse, ')'),
		label = str_replace(label, '\\)', mouse)) %>% 
	group_by(label) %>% 
	mutate(median_p = median(Cleared)) %>% 
	filter(0.45 < median_p, median_p < 0.55) %>% 
	select(Group, sample_label = label, median_p) %>% 
	unique
# plot relative abundance of top otus for potentially misclassified features
rel_abund_misclass_samples <- top_otu_abundance %>% 
	inner_join(summarise(data_full, new_auc = median(new_auc)), by = c('label' = 'names')) %>% 
	right_join(misclass_samples, by = c('Group')) %>% 
	mutate(sample_label = gsub(' \\(', '\\\n\\(', sample_label)) %>% 
	ggplot(aes(x = fct_reorder(label, -new_auc), y = abundance, color = clearance)) + 
		geom_point() + 
		scale_y_log10() +
		coord_flip() + 
		facet_wrap(.~reorder(sample_label, median_p), nrow = 1) + 
		theme_bw() +
		labs(y = 'Percent Relative Abundance', x = NULL, color = NULL) + 
		theme(legend.position = 'none',
			panel.grid.minor.y = element_blank(),
			axis.text.y = element_markdown())

# -------------------------------------------------------------------->


######################################################################
#-------------------Plot facultative anaerobes---------------------- #
######################################################################

mice <- l2_sample_perf_df %>% 
	left_join(full_meta_df, by = c('Group' = 'group')) %>% 
	pull(mouse_id) %>% 
	unique

fac_anaerobe_plot <- full_meta_df %>% 
	filter(mouse_id %in% mice) %>% 
	select(mouse_id, day, clearance, abx, CFU, Group = group) %>% 
	inner_join(select(shared_df, Group, Otu000010, Otu000011), by = 'Group') %>% 
	pivot_longer(cols = starts_with('Otu0'), names_to = 'otu', values_to = 'abundance') %>% 
	inner_join(label_df, by = c('otu' = 'key')) %>% 
	mutate(abundance = (100 * abundance/total_abundance) + .05) %>% 
	filter(day > 0) %>% 
	ggplot(aes(x = day, y = abundance, color = label)) + 
		stat_summary(fun = function(x) median(x), geom = "line", size = 3) + # Median darker
		geom_line(aes(group = interaction(mouse_id, label)), alpha = 0.3) + 
		scale_x_continuous(breaks = 1:10) + 
		scale_y_log10() + 
		facet_grid(abx~clearance) + 
		theme_bw() + 
		theme(legend.position = c(0.8, 0.5),
			panel.grid.minor = element_blank()) +
		labs(x = 'Day', y = 'Relative Abundance', color = NULL)

cfu_plot <- full_meta_df %>% 
	filter(mouse_id %in% mice) %>% 
	select(mouse_id, day, clearance, abx, CFU, Group = group) %>% 
	filter(day > 0) %>% 
	ggplot(aes(x = day, y = CFU)) + 
		stat_summary(fun = function(x) median(x), geom = "line", size = 3) + # Median darker
		geom_line(aes(group = mouse_id), alpha = 0.3) + 
		scale_x_continuous(breaks = 1:10) + 
		scale_y_log10() + 
		facet_grid(abx~clearance) + 
		theme_bw() + 
		theme(legend.position = c(0.8, 0.5),
			panel.grid.minor = element_blank()) + 
		labs(x = 'Day', y = expression(italic('C. difficile')~' CFU'))


# -------------------------------------------------------------------->



######################################################################
#-----------------------Save figure -------------------------------- #
######################################################################
#combine with cowplot


ggsave(paste0("results/figures/Figure_4.jpg"), 
	plot = plot_grid(
			plot_grid(NULL, perm_imp_plot, NULL, rel_heights = c(1,45,2), ncol = 1),
			top_otu_abundance_plot,
			labels = c('A', 'B'), nrow = 1),
	width = 12, height = 10, units="in")

ggsave(paste0("results/figures/Figure_S4.jpg"), 
	plot = model_perf_plot,
	width = 6, height = 6, units="in")

ggsave(paste0("results/figures/Figure_S4_L2_Logistic_Regression_sample_dist_", model_dir, ".jpg"), 
	plot = plot_grid(plot_grid(NULL, perf_by_sample_plot, rel_heights = c(1, 19), ncol = 1), 
			rel_abund_misclass_samples, 
		labels = c('A', 'B'), rel_widths = c(1,2)), 
	width = 16, height = 8, units="in")

ggsave(paste0("results/figures/Figure_S4_fac_anaerobes_", model_dir, ".jpg"), 
	plot = plot_grid(fac_anaerobe_plot, cfu_plot, 
		labels = c('A', 'B')),  
	width = 16, height = 8, units="in")
