library(tidyverse)

cost_5050_l2 <- read_csv('data/process/otu_5050/combined_all_hp_results_L2_Logistic_Regression.csv')
cost_5050_rf <- read_csv('data/process/otu_5050/combined_all_hp_results_Random_Forest.csv')
data_df <- bind_rows(
    select(mutate(cost_5050_l2, hp = 'cost', params = cost), 
    	hp, params, ROC),
	select(mutate(cost_5050_rf, hp = 'mtry', params = mtry), 
		hp, params, ROC))

plot <- data_df %>% 
     ggplot(aes(x = params, y = ROC, group = params, color = hp)) + 
         geom_boxplot() + 
         #scale_x_log10() + 
         facet_grid(.~hp, scales = 'free')

cost_5050_rf %>% 
     ggplot(aes(x = mtry, y = ROC, group = mtry)) + 
         geom_boxplot() 


ggsave('scratch/fig_S4_tuning_C.jpg', plot)