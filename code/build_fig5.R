##############
#
# run script to generate plots for Figure 4
#	How do communities predict/classify clearance/colonization?
# 
# Nick Lesniak 04-28-2020
#
#  need files:
#	data/process/abx_cdiff_taxonomy_clean.tsv
#	data/process/abx_cdiff_metadata_clean.txt
#	data/process/l2_otu/combined_best_hp_results_L2_Logistic_Regression.csv
#	data/process/l2_otu/combined_all_sample_results_L2_Logistic_Regression.csv
#	data/process/l2_otu/combined_L2_Logistic_Regression_feature_ranking.tsv
#	data/process/l2_otu/combined_best_hp_results_L2_Logistic_Regression.csv
#	data/process/l2_otu/L2_Logistic_Regression_non_cor_importance.tsv
#	data/process/l2_otu/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv
#	data/mothur/sample.final.0.03.subsample.shared
#	code/R/functions.R
#
##############

library(tidyverse)
library(cowplot)
library(pROC)
library(ggtext)

source("code/R/functions.R")

#args <- commandArgs(trailingOnly = TRUE)
model_dir <- 'l2_otu'

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

# colors for each antibiotic
abx_color <- tibble(abx = c('Clindamycin', 'Cefoperazone', 'Streptomycin'),
	color = c('#A40019', '#3A9CBC', '#D37A1F'))


#wilcox.test(l2_performance_df$cv_aucs,l2_performance_df$test_aucs)$p.value
# plot the aucs
model_perf_plot <- l2_performance_df %>% 
	rename(`Test\nAUC` = test_aucs, `CV\nAUC` = cv_aucs) %>% 
	pivot_longer(cols = contains('AUC'), values_to = 'AUC', names_to = 'validation', ) %>% 
	ggplot(aes(x = validation, y = AUC)) +
		geom_boxplot(alpha=0.5, fatten = 2) +
		geom_hline(yintercept = 0.5, linetype="dashed") +
		coord_cartesian(ylim = c(0.4, 1)) +
		theme_bw() +
		theme(legend.position='none',
		      panel.grid.major.y = element_blank(),
		      panel.grid.minor = element_blank(),
		      panel.background = element_blank()) + 
		labs(x = NULL, y = "AUROC") 
		


######################################################################
#--------------Run the functions and plot feature ranks ----------#
######################################################################

label_df <- bind_rows(tax_df, select(meta_df, -mouse, -Group), abx_df) %>% 
	mutate(label = paste0('*', label),
		label = gsub(' \\(', '* \\(', label))

# -------------------------------------------------------------------->


######################################################################
#------------- Plot L2 Log Reg Permutation Importance -------------- #
######################################################################

# Read in the cvAUCs, testAUCs for 100 splits as base test_aucs
logit <- read_files(paste0("data/process/", model_dir,
	"/combined_best_hp_results_L2_Logistic_Regression.csv"))
# ----------------------------------------------------------------------------->


# --------  Get the top OTUs that have the largest impact on AUROC ---------->

# Define the function to get the  most important top OTUs
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
perm_imp_plot <- ggplot(data_full, aes(fct_reorder(names, -new_auc), new_auc)) +
	annotate('ribbon', x = c(-Inf, Inf), ymin= lowerq, ymax= upperq, fill = 'green4', alpha = 0.2) +
	geom_hline(yintercept = data_base_medians$imp , color = 'green4') +
	geom_boxplot(fill = 'white', width = 0.5) +
	coord_flip() +
	theme_bw() +
	labs(y = "AUROC", x = NULL) + 
		#x = expression(paste(L[2], "-regularized logistic regression"))) +
	theme(legend.position="none",
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.text.y=element_markdown())

# -------------------------------------------------------------------->

######################################################################
#------------------  Plot feature coefficients  -------------------- #
######################################################################
otu_order <- data_full %>% 
	group_by(names) %>% 
	summarise(new_auc = median(new_auc))

logit_imp <- read_tsv(paste0("data/process/", model_dir, 
	"/combined_L2_Logistic_Regression_feature_ranking.tsv")) %>% 
	left_join(label_df, by = c('key')) %>% 
	group_by(key) %>% 
	mutate(sign = ifelse(names(which.max(table(sign))) == 'negative', 'Cleared', 'Colonized'),
		p_cleared =  1 - (exp(value)/(1+exp(value))), # convert logit of remaining colonized to probability of clearance 
		OR_cleared = p_cleared/(1-p_cleared),
		logit_cleared = log(OR_cleared)) %>% # convert probability to odds ratio
	inner_join(otu_order, by = c('label' = 'names'))  %>% 
	ungroup

axis_colors <- logit_imp %>% 
	select(label, sign, new_auc) %>% 
	distinct() %>% 
	mutate(sign = ifelse(sign == 'Cleared', '#C14642', '#008E94'),
		label = fct_reorder(label, -new_auc)) %>% 
	arrange(label) %>% 
	pull(sign)

coef_plot <- logit_imp %>% 
	ggplot(aes(x = fct_reorder(label, -new_auc), y = OR_cleared, color = sign)) +
      geom_boxplot(width = 0.5, show.legend = F) +
      geom_point(size = NA, shape = 15) +
      geom_hline(yintercept = 1, linetype='dashed') + 
      coord_flip() +
      theme_bw() +
      theme(legend.title = element_blank(),
          legend.position = c(0.85,0.05),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) + 
	  guides(colour = guide_legend(override.aes = list(size = 5))) + 
      labs(y = 'Odds ratio', x = NULL) 
# -------------------------------------------------------------------->

######################################################################
#------------------ Do the features make sense? -------------------- #
######################################################################
# plot abundance by outcome
shared_file <- 'data/mothur/sample.final.0.03.subsample.shared'
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
plot_abundance <- function(antibiotic){
	abx_col <- pull(filter(abx_color, abx == antibiotic), color)

	p <- top_otu_abundance %>% 
		filter(abx == antibiotic) %>% 
		mutate(abx = factor(abx, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin'))) %>% 
		inner_join(otu_order, by = c('label' = 'names')) %>% # add perm importance aucs to order otus
		ggplot(aes(x = fct_reorder(label, -new_auc), y = abundance, color = clearance)) + 
			geom_point(position = position_jitterdodge(jitter.width = 0.25), alpha = .3) + 
			scale_y_log10(limits = c(0.04,100),
			   		breaks = c(0.01, 0.1, 1, 10, 100),
			   		labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) + 
			coord_flip() + 
			facet_wrap(.~abx) + 
			theme_bw() +
			labs(x = NULL, color = NULL) + 
			theme(panel.grid.minor = element_blank(),
				axis.text.x = element_markdown(),
				strip.background = element_rect(fill = abx_col),
				strip.text = element_text(color = 'white'))
	if(antibiotic == 'Clindamycin'){
		p + theme(legend.position = 'none',
		          axis.ticks.y = element_blank(),
		          axis.text.y = element_blank()) + 
			labs(y = NULL)
	} else if(antibiotic == 'Cefoperazone') {
		p + theme(legend.position = 'none',
				axis.ticks.y = element_blank(),
				axis.text.y = element_blank()) + 
			labs(y = 'Relative abundance (%)') +
			guides(colour = guide_legend(override.aes = list(alpha = 1)))
	} else {
		p + theme(legend.position = 'none',
				axis.ticks.y = element_blank(),
				axis.text.y = element_blank()) + 
			labs(y = NULL)
	}
}

top_otu_abundance_plot <- plot_grid(
	plot_grid(plot_abundance('Clindamycin'), NULL, ncol = 1, rel_heights = c(50, 1)), 
	plot_abundance('Cefoperazone'), 
	plot_grid(plot_abundance('Streptomycin'), NULL, ncol = 1, rel_heights = c(50, 1)), 
	nrow = 1)

######################################################################
#-----------------------Save figure -------------------------------- #
######################################################################
#combine with cowplot


ggsave(paste0("submission/figure_5.tiff"), 
	plot = plot_grid(
	  plot_grid(NULL, NULL, NULL, labels = c('A', 'B', 'C'), nrow = 1, rel_widths = c(5,4,4)),
	  plot_grid(
			plot_grid(NULL, perm_imp_plot, rel_heights = c(1,45), ncol = 1),
			plot_grid(NULL, coef_plot, rel_heights = c(1,45), ncol = 1),
			top_otu_abundance_plot,
			nrow = 1, rel_widths = c(5,4,4)),
	  ncol = 1, rel_heights = c(1, 40)),
	width = 18, height = 10, units="in", compression = 'lzw')

ggsave(paste0("submission/figure_S5.tiff"), 
	plot = model_perf_plot,
	width = 6, height = 6, units="in", compression = 'lzw')