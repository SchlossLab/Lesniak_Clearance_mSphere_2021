
library(tidyr)
library(dplyr)
#library(ggplot2)
#library(scatterplot3d)
library(plotly)
library(rEDM)
#library(patchwork)

lagged <- 1 # time offset

input_file <- read.table('data/process/abx_cdiff_metadata_clean.txt',
	header = T, stringsAsFactors = F)

# create lagged dataframe
input_df <- input_file %>%
	select(CFU, cage, mouse, abx, day) %>%  
	filter(abx == 'cef', cage == 600, day > 0) %>%
	arrange(mouse, day) %>%
	group_by(mouse) %>%
	mutate(CFU_0lag = log10(CFU + 1),
		CFU_1lag = c(NA, CFU_0lag[-max(day)]),
		CFU_2lag = c(NA, CFU_1lag[-max(day)]),
		CFU_3lag = c(NA, CFU_2lag[-max(day)]))

plot_ly(test_df, x = ~CFU_0lag, y = ~CFU_1lag, z=~CFU_2lag,
	mode = 'lines', type = 'scatter3d', color = ~mouse,
	line = list(width = 6)) 

cbind(select(filter(test_df, mouse == 1), x1 = CFU_0lag, y1 = CFU_1lag, z1 = CFU_2lag, day),
	select(filter(test_df, mouse == 2), x2 = CFU_0lag, y2 = CFU_1lag, z2 = CFU_2lag),
	select(filter(test_df, mouse == 3), x3 = CFU_0lag, y3 = CFU_1lag, z3 = CFU_2lag)
) %>%
	plot_ly(x = ~x1, y = ~y1, z = ~z1, type = 'scatter3d', mode = 'lines',
		colors = 'Set1',
        line = list(color = 'orange', width = 6)) %>%
	add_trace(x = ~x2, y = ~y2, z = ~z2,
		colors = 'Set2',
		line = list(color = 'blue', width = 6)) %>%
	add_trace(x = ~x3, y = ~y3, z = ~z3,
		line = list(color = 'red', width = 6))

#	#community dynamics plot
#	input_df %>% 
#		gather(species, abundance, -epoch) %>% 
#		filter(species %in% c('Bulbasaur', 'Charmander', 'Squirtle')) %>% 
#		ggplot(aes(x = epoch, y = abundance, color = species)) +
#			geom_line() +
#			ylim(0,NA)
#	# attractor manifold (M)
#	attractor_manifold_int <- plot_ly(input_df, 
#		x = ~Bulbasaur, y = ~Charmander, z = ~Squirtle, 
#		mode = 'lines', type = 'scatter3d')
#	png('exploratory/scratch/temp_plot.png', width = 900, height = 900)
#	with(input_df, scatterplot3d(Bulbasaur, Charmander, Squirtle, 
#		type = 'l', angle = 150))
#	dev.off()
#	g <- grid::rasterGrob(png::readPNG('exploratory/scratch/temp_plot.png'), interpolate = T)
#	gg_attractor_manifold <- ggplot(input_df, aes(x = Bulbasaur, y = Charmander)) + 
#		annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
#		theme(axis.title = element_blank())
#	# plots to test for non-linear dynamical signature
#	# plot manifold
	temporal_plot <- test_df %>% 
		ggplot(aes(x = day, y = CFU, color = as.factor(mouse))) +
			geom_line() +
			ylim(0,NA) + scale_y_log10() +
			theme_classic() + 
			theme(legend.position = c(0.9, 0.9), 
				legend.background = element_rect(size=0.25, linetype="solid",
					colour ="black")) + 
			labs(title = paste(unique(test_df$abx)), color = 'Mouse', x = 'Day', y = 'CFU (log10)')
#	shadow_manifold_plot <- plot_ly(input_df, 
#		x = ~Charmander, y = ~Charmander_1lag, z = ~Charmander_2lag,
#		mode = 'lines', type = 'scatter3d')
#	
#	#simplex_plots - identify best dimension
#	set.seed(041618)

simplex_c_plot <- c()
for(mouse_number in (unique(test_df$mouse))){
	for(portion in seq(0.1, 0.9, 0.1)){
	input_df <- filter(test_df, mouse == mouse_number) %>%	ungroup() %>%
		select(day, CFU) %>% arrange(day) %>%
		data.frame
	train <- c(1, round(portion* nrow(input_df))) # time points to subset for training
	test <- c( ( round(portion* nrow(input_df)) + 1 ), nrow(input_df)) # time points to subset for testing
	
	simplex_out <- #do.call(rbind, lapply('CFU', function(var) {
#	    cbind(
	    	simplex(input_df, E = 1:10, lib = train, pred = test)#,
#	    	species = var)
#	}))
	simplex_c_plot <- rbind(simplex_c_plot,
		cbind(select(simplex_out, E, rho), mouse = mouse_number, train_portion = portion, sample = paste(mouse_number, portion, sep = '_')))
}}
#embedding_dim_plot <- 
simplex_c_plot %>% 
	ggplot(aes(x = E, y = rho, group = sample)) +
		geom_point() + geom_line() + 
		facet_grid(mouse ~ train_portion) +
		labs(x = "Embedding Dimension (E)", y = "Forecast Skill (rho)",
			title = 'Simplex') +
		theme_bw()
#	
#best_E <- 
simplex_c_plot %>% 
	filter(!is.na(rho)) %>%
	group_by(mouse) %>% 
	mutate(best = max(rho)) %>% 
	filter(best == rho, rho > 0) %>% ungroup #%>% 
	#with(., setNames(E, as.character(mouse)))

best_E <- 3
train <- c(1, round(0.5 * max(test_df$day))) # time points to subset for training
test <- c( round(0.5 * max(test_df$day) + 1 ), max(test_df$day)) # time points to subset for testing
#	
#	#s-maps - identify non-linearity
#	smap_out <- do.call(rbind, lapply(names(input_df)[2:4], function(var) {
#	    cbind(
smap_plot <- c()
for(mouse_number in (unique(test_df$mouse))){
	input_df <- filter(test_df, mouse == mouse_number) %>%	ungroup() %>%
		select(day, CFU) %>% arrange(day) %>%
		data.frame
	smap_out <- s_map(input_df, E = best_E, lib = train, pred = test)
	smap_plot <- rbind(smap_plot,
		cbind(select(smap_out, theta, rho), mouse = mouse_number))
}
#	    	species = var)
#	}))
#	
#	nonlinear_plot <- 
smap_plot %>% 
	ggplot(aes(x = theta, y = rho)) +
		geom_line() + facet_grid(mouse~., scales = 'free_y') +
		labs(x = "Nonlinearity (theta)", y = "Forecast Skill (rho)",
			title = "S-Map") +
		theme_bw()
#	
#	#forecast using simplex/s-map
#	
#	pred_obs_plot <- c()
#	for(var in species_list){
#		target <- var
#		cols <- names(input_df)[startsWith(names(input_df), var)]
#		block_lnlp_output <- block_lnlp(input_df, lib = train, pred = test, columns = cols, 
#	    	target_column = target, stats_only = FALSE, first_column_time = TRUE)
#		pred_obs_plot <- rbind(pred_obs_plot, data.frame(
#			observed =  block_lnlp_output$model_output[[1]]$obs,
#			predicted = block_lnlp_output$model_output[[1]]$pred,
#			rho = block_lnlp_output$rho,
#			species = target, input_variables = paste(cols, collapse = '__'))) 
#	}
#	
#	forecast_plot <- pred_obs_plot %>% 
#			ggplot(aes(x = observed, y = predicted, color = as.factor(rho))) + 
#			geom_point(alpha = 0.1) + facet_grid(species~., scales = 'free_y') +
#			geom_abline(slope = 1, intercept = 0, linetype = 2) +
#			geom_smooth(method = 'lm') +  
#			labs(x = "Observed", y = "Predicted", color = 'rho') +
#			theme_bw() + theme(legend.position = c(0.9, 0.2))
#	
#	# recover influence via cross mapping w/simplex
#	library_sizes <- c(seq(5, 55, by = 2), seq(55, 400, by = 50))
#	
#	xmap_plot <- c()
#	for(var in species_list){
#		lib_col <- var
#		for(tar in species_list[species_list != var]){
#			target_col <- tar
#			library_xmap_target <- ccm(input_df, lib = train, pred = test, 
#				lib_column = lib_col, target_column = target_col, 
#				E = best_E[lib_col], lib_sizes = library_sizes, silent = TRUE)
#			xmap_means <- ccm_means(library_xmap_target)
#			xmap_plot <- rbind(xmap_plot, data.frame(
#				library_size = xmap_means$lib_size, 
#				rho = pmax(0, xmap_means$rho),
#				comparison = paste0(lib_col, '_xmap_', target_col)))
#		}
#	}
#	
#	crossmap_plot <- xmap_plot %>% 
#		ggplot(aes(x = library_size, y = rho, color = comparison)) + 
#			geom_line() + 
#			labs(x = "Library Size", y = "Cross Map Skill (rho)") +
#			theme_bw() +
#			theme(legend.position = c(0.75, 0.2), legend.title = element_blank(),
#				legend.background = element_rect(size=0.25, linetype="solid",
#					colour ="black"))
#	
#	#with(delta_df[between(delta_df$epoch, 10, 310), ],
#	#	lines3d(Charmander, char_1, char_2, 
#	#		xlab = 'Bulb', ylab = 'Char', zlab = 'Sqrtl',
#	#		width = 1000, size=3, col = 'red'))
#	#
#	#with(delta_df[between(delta_df$epoch, 10, 310), ],
#	#	lines3d(Squirtle, sqtl_1, sqtl_2, 
#	#		xlab = 'Bulb', ylab = 'Char', zlab = 'Sqrtl',
#	#		width = 1000, size=3, col = 'blue'))
#	
#	
#	ggsave('exploratory/scratch/CCM_rps.jpg',
#		( ( gg_attractor_manifold / temporal_plot ) | ( embedding_dim_plot + nonlinear_plot) )  /
#		(forecast_plot  + crossmap_plot),
#		height = 14, width = 12)#	