library(cowplot)
library(rEDM)
library(tidyverse)

input_values <- commandArgs(TRUE)
run_set <- as.numeric(input_values[1])
print(paste0('Running set ', run_set))
input_file <- as.character(input_values[2])
#input_file <- 'data/process/ccm_otu_data.txt'
#input_file <- 'data/process/ccm_validation_data.txt'
#input_file <- 'data/process/bucci/ccm_bucci_data.txt'
otu_df   <- read.table(input_file, header = T, stringsAsFactors = F) %>% 
  mutate(otu_feature = gsub('_first', '', otu_feature)) %>%
  filter(differenced == 'first')
iters <- 500
seed <- 062818
treatment_subset <- unique(otu_df$treatment)[run_set]
print(paste0('Running set ', run_set, ' - Treatment ', treatment_subset))
otu_df <- filter(otu_df, treatment %in% treatment_subset)
save_dir <- paste0('data/process/ccm/', treatment_subset, '/interactions')
ifelse(!dir.exists(save_dir), 
  dir.create(save_dir), 
  print(paste0(save_dir, ' directory ready')))

sample_list <- unique(otu_df$unique_id)
otu_list <- unique(otu_df$taxa)
## import results from individual taxa test
## check for embedding
#load_files <- c('simplex_embedding_first_differenced.txt')
#for(i in load_files){
#  if(!file.exists(paste0(save_dir, '/../', i))){ 
#    stop(paste('Run previous script or locate', i))
#    }
#  }
##  import embedding
#embedding <- read.table(paste0(save_dir, '/../', load_files[1]), 
#  header = T, stringsAsFactors = F) %>% 
#    mutate(taxa = gsub('_first', '', taxa))

# if validation, import known interactions
true_interactions <- read.table(paste0('data/process/validation/validation_interaction_matrix_', 
  treatment_subset, '.txt')) %>% 
  mutate(affected_otu = otu_list) %>% 
  gather(affector_otu, actual_strength, -affected_otu) %>% 
  mutate(affector_otu =  gsub('V', 'OTU_', affector_otu))

run_smap <- function(i, j){
  composite_ts <- otu_df %>% 
    filter(taxa %in% c(i, j)) %>% 
    select(taxa, day, unique_id, normalized_abundance) %>% 
    spread(taxa, normalized_abundance) %>% 
    arrange(unique_id, day)   
  
  data_by_plot <- split(composite_ts, composite_ts$unique_id)
  segments_end <- cumsum(sapply(data_by_plot, NROW))
  segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
  segments <- cbind(segments_begin, segments_end) 
  pred_num <- round(0.2 * length(sample_list))
  lib_segments <- segments[-tail(1:NROW(segments), pred_num), ]
  pred_segments <- segments[tail(1:NROW(segments), pred_num), ]

  univariate_ts <- bind_cols(composite_ts[, c('day', 'unique_id')], 
    make_block(composite_ts[,i],
      lib = segments))#,

  # check mae for best theta
  test_theta <- list()
  for(iter in 1:iters){
    rndpred <- sample(1:length(sample_list), pred_num)
    pred_order <- data.frame(
      unique_id = c(sample(sample_list[-rndpred], replace = T), sample_list[rndpred]),
      order = 1:NROW(sample_list), 
      stringsAsFactors = F)
    rnd_composite_ts <- right_join(composite_ts, pred_order, by = 'unique_id')
    rnd_univariate_ts <- right_join(univariate_ts, pred_order, by = 'unique_id')

    uni_theta <- block_lnlp(select(univariate_ts, contains('col1')),
      method = 's-map',
      num_neighbors = 0,
      lib = lib_segments, pred = pred_segments,
      theta = c(0, 1e-04, 3e-04, 0.001,
        0.003, 0.01, 0.03, 0.1,
                0.3, 0.5, 0.75, 1, 1.5,
                2, 3, 4, 6, 8),
      target_column = 'col1',
      silent = T)
    theta_run <- block_lnlp(select(composite_ts, one_of(i, j)),
      method = 's-map',
      num_neighbors = 0,
      lib = lib_segments, pred = pred_segments,
      theta = c(0, 1e-04, 3e-04, 0.001,
        0.003, 0.01, 0.03, 0.1,
                0.3, 0.5, 0.75, 1, 1.5,
                2, 3, 4, 6, 8),
      target_column = i,
      silent = T)
    test_theta[[iter]] <- bind_rows(
      data.frame(theta_run, model_type = 'multivariate', stringsAsFactors = F),
      data.frame(uni_theta, model_type = 'univariate', stringsAsFactors = F))
  }

  best_theta <- do.call('rbind', test_theta) %>% 
    group_by(model_type, theta) %>% 
    summarise(median_mae = median(mae, na.rm = T)) %>% 
    filter(median_mae == min(median_mae)) 

  # run smap with best theta
  # partial derivaties of multivariate smap approximates interactions
  interaction_smap <- list()
  for(iter in 1:iters){
    pred_order <- data.frame(unique_id = sample(sample_list, replace = T),
      order = 1:NROW(sample_list), stringsAsFactors = F)
    rnd_composite_ts <- right_join(composite_ts, pred_order, by = 'unique_id')
    rnd_univariate_ts <- right_join(univariate_ts, pred_order, by = 'unique_id')
    smap_multi <-  block_lnlp(select(rnd_composite_ts, one_of(i, j)),
      lib = segments, pred = segments,
      method = "s-map",
      num_neighbors = 0, 
      theta = pull(filter(best_theta, model_type == 'multivariate'), theta),
      target_column = i,
      silent = T,
      save_smap_coefficients = T) # save S-map coefficients
    smap_uni <-  block_lnlp(select(univariate_ts, contains('col1')),
      lib = segments, pred = segments,
      method = "s-map",
      num_neighbors = 0, 
      theta = pull(filter(best_theta, model_type == 'univariate'), theta),
      target_column = 'col1',
      silent = T,
      save_smap_coefficients = T) # save S-map coefficients
    smap_coef_df <- smap_multi$smap_coefficients[[1]]
    colnames(smap_coef_df) <- c(paste0('d', i, '_d', i), paste0('d', i, '_d', j), 'intercept')
    interaction_smap[[iter]] <- bind_rows(
      data.frame(bind_cols(smap_multi$model_output[[1]], smap_coef_df,
          right_join(composite_ts, pred_order, by = 'unique_id')), 
        embed = 'multi', run = iter, mae = smap_multi$mae, stringsAsFactors = F),
      data.frame(smap_uni$model_output[[1]], embed = 'uni', run = iter, 
        mae = smap_uni$mae, stringsAsFactors = F))    
  }
  output <- list(model_output = do.call('rbind', interaction_smap))
  output[['mae']] <- cbind(data.frame(affected_otu = i, affector_otu = j,
        mae = output[['model_output']] %>%
          filter(embed == 'multi') %>%
          pull(mae) %>%
          median, stringsAsFactors = F),
      output[['model_output']] %>% 
        filter(embed == 'multi') %>%
        gather(interaction, strength, contains('dOTU')) %>%
        group_by(interaction) %>% 
        summarise(median_interaction = median(strength, na.rm = T),
          lower_quartile = quantile(strength, na.rm = T)['25%'],
          upper_quartile = quantile(strength, na.rm = T)['75%'])
        )
  return(output)
}

# for each otu
for(i in otu_list){
  smap_output <- lapply(otu_list[otu_list != i], function(j) run_smap(i, j))
  mae_list <- c()
  for(n in 1:length(otu_list[otu_list != i])){
    mae_list <- rbind(mae_list, smap_output[[n]]$mae)
  }

# for all other otus, test for prediction by s-map
  step_perf <- list()
  step_mae <- c()

step_mae
for( in 
step_increase[[]] %>%
  filter(embed == 'multi') %>%
  pull(mae) %>%
  median
  table

tmp <- do.call('rbind', step_increase)
tmp %>%
    gather(interaction, interaction_strength, contains('dOTU')) %>%
    gather(OTU_time, OTU_abundance, contains('OTU'))
  filter(mae == min(mae)) %>%
  pull(run) %>%
  unique

  write.table(interaction_smap, paste0(save_dir, '/interactions_w_', i, '.txt'), 
  quote = F, row.names = F)

  # comparison to univariate is limited to late time points 
  # so  limited interpretation
  #interaction_smap %>% 
  # ggplot(aes(x = obs, y = pred, color = embed)) + 
  #   geom_point(alpha = .2) + 
  #   theme_bw() + 
  #     coord_cartesian(ylim=c(-2, 2), xlim=c(-2, 2))
  #interaction_smap %>% 
  # group_by(embed) %>% 
  # summarise(mean = mean(obs, na.rm = T), median = median(obs, na.rm = T),
  #   pval = cor.test(.$obs, .$pred, method = 'spearman')$p.val,
  #   cor = cor.test(.$obs, .$pred, method = 'spearman')$estimate)

  ## Time series of fluctuating interaction strength
  order <- expand.grid(unique_id = sample_list, day = 0:10)  %>% 
    arrange(unique_id, day) %>% 
    mutate(epochs = 1:NROW(.)) #%>% 
    #right_join(expand.grid(epochs = 1:NROW(.), run = 1:500)) %>% 
    #arrange(unique_id, run, day)
  # Plot all partial derivatives
  interaction_temporal_plot <- interaction_smap %>% 
    filter(embed == 'multi') %>% 
    right_join(order, by = c('day', 'unique_id')) %>% 
    gather(interaction, strength, one_of(paste0('d', i, '_d', j))) %>% 
    mutate(interaction = gsub('tu0*', 'TU', interaction)) %>% 
    ggplot(aes(x = epochs, y = strength)) + 
      stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
        alpha = 0.2, fun.args =(conf.int = 0.5), aes(fill = interaction)) + 
      stat_summary(aes(color = interaction), fun.y = median, geom = 'line') + 
      theme_bw() +
      theme(legend.position="top", legend.title=element_blank()) +
      labs(title = paste0('Strength of effect of OTUs on ', gsub('tu0*', 'TU', i)),
        subtitle = 'Median with IQR shown')

  interaction_dist_plot <- interaction_smap %>% 
    filter(embed == 'multi') %>% 
    right_join(order, by = c('day', 'unique_id')) %>% 
    gather(interaction, strength, one_of(paste0('d', i, '_d', j))) %>% 
    mutate(interaction = gsub('^.*_d', '', interaction)) %>% 
    inner_join(filter(true_interactions, affected_otu == i), ) %>% 
    ggplot(aes(x = interaction, y = strength, fill = interaction)) + 
      geom_hline(yintercept = 0, linetype = 'dashed') + 
      geom_violin(alpha = 0.5) + 
      geom_boxplot(aes(y = actual_strength), color = 'red') + 
      theme_bw() +
      theme(legend.position='none') +
      labs(title = paste0('Strength of effect of OTUs on ', gsub('tu0*', 'TU', i)),
        subtitle = 'Median with IQR shown')
  
# dynamics_plot <- shared_file %>% 
#     select(cage, mouse, day, one_of(i, j)) %>% 
#     mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
#     right_join(order) %>% 
#     #mutate(time = 1:nrow(.)) %>% 
#     gather(bacteria, counts, one_of(i, j)) %>% 
#       ggplot(aes(x = epochs, y = counts, color = as.factor(cage), group = interaction(cage, mouse))) + 
#         geom_line(alpha = 0.4) + 
#         geom_point() + 
#         facet_grid(bacteria~., scales = 'free_y') +
#         theme_bw() + 
#         labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
#           subtitle = 'Temporal Dynamics - Colored by mouse', 
#           title = paste(i, 'is driven by', paste(j, collapse = ', '))) + 
#         theme_bw(base_size = 8) + 
#         theme(legend.position = 'none', panel.grid.minor = element_blank())

  ggsave(filename = paste0(save_dir, '/interactions_w_', i, '.jpg'), 
    plot = plot_grid(interaction_temporal_plot, interaction_dist_plot, nrow = 2, align = 'v'),
    width = 14, height = 14, device = 'jpeg')
  

#}
}
# output file with taxa, interaction direction, interaction strength



#GenerateDdbDtg[data_] := Module[{dbLength, xs, ys, timeStamps,
#   fordb, dbtk, rndStp, xdb, ydb, fortg, xtg, ytg},
#  {xs, ys} = data;
#  timeStamps = Range[Length[xs]];
#  dbLength = IntegerPart[Length[timeStamps]*0.5] + 1;
#  dbtk = Min[dbLength, 100];
#  rndStp = RandomSample[timeStamps];
#  fordb = Take[rndStp, dbtk];
#  xdb = xs[[fordb]];
#  ydb = ys[[fordb]];
#  fortg = Drop[rndStp, dbtk];
#  xtg = xs[[fortg]];
#  ytg = ys[[fortg]];
#  {xdb, ydb, xtg, ytg}]
#
#Smap <- function(xdb, ydb, xtg, ytg, theta) := 
# Module[{yStar, ck, weights, as, bs, cs, sxdb, norms, ds},
#  {yStar, ck} = Transpose[(
#       sxdb = Transpose[Transpose[xdb] - INPUT_];
#       norms = Norm /@ sxdb;
#       ds = Mean[norms];
#       weights = N[Exp[-theta*(norms/ds)]];
#       as = weights*xdb;
#       bs = weights*ydb;
#       cs = PseudoInverse[as].bs;
#       {cs.INPUT_, cs}) & /@ xtg];
#  {Total[(yStar - ytg)^2]/Length[ytg],
#   Mean[ck](*ctilde*)}]
#
#(*Smap result with the best theta*)
#
#Bestsmap[iactive_, i_, data_] := 
# Module[{activePositions, xdb, ydb, xtg, ytg,
#   bestResult, mse, cTilde, theta, itest, sres, best, iscale, step, 
#   neighb, newtest,
#   newres},
#  activePositions = Flatten[Position[iactive, 1]];
#  {xdb, ydb, xtg, ytg} = {data[[1, All, activePositions]], 
#    data[[2, All, i]],
#    data[[3, All, activePositions]], data[[4, All, i]]};
#  (*minimization of theta*)
#  itest = Range[0, 8];
#  sres = Sort[
#    Table[{Smap[xdb, ydb, xtg, ytg, theta], theta}, {theta, itest}]];
#  best = sres[[1]];
#  iscale = 1; step = 0;
#  While[step < 
#     5 && (step != 2 || best[[2]] >= 0.5) && (step != 1 || 
#      best[[2]] >= 0.5),
#   (*0.5以下の場合は無視する設定.これが特にノイズが有る場合などに,
#   影響しないことは確認するべき*)
#   step++;
#   neighb = Select[sres, (Abs[INPUT_[[2]] - best[[2]]] == iscale &)];
#   newtest = (best[[2]] + neighb[[All, 2]])/2;
#   newres = 
#    Sort[Table[{Smap[xdb, ydb, xtg, ytg, theta], theta}, {theta, 
#       newtest}]];
#   sres = Sort[Join[{best}, neighb, newres]];
#   best = sres[[1]];
#   iscale = iscale/2
#   ];
#  bestResult = best;
#  (*MSE and Overscript[c, ~]i (ci tilde)*)
#  
#  theta = N[bestResult[[2]]];
#  mse = bestResult[[1, 1]];
#  cTilde = bestResult[[1, 2]];
#  {mse,
#   ReplacePart[iactive, 
#    MapThread[(INPUT_1 -> INPUT_2) &, {activePositions, cTilde}]](*cij^**), 
#   theta}
#  ]
#
#SSM[timeseries_, i_, qc_, baggingIteration_] := (
#  baggingResult = ParallelTable[
#    numberofSpecies = Length[timeseries[[1]]];
#    ddbDtg = GenerateDdbDtg[GenerateXY[timeseries]];
#    iactive = 
#     Append[ReplacePart[Table[0, {numberofSpecies}], i -> 1], 1];
#    initialResult = Bestsmap[iactive, i, ddbDtg];
#    msePrev = initialResult[[1]];
#    ciStar = initialResult[[2]];
#    thetaBest = initialResult[[3]];
#    links = 1; q = 99;
#    While[links < numberofSpecies - 1 && q > qc,
#     newIactives = 
#      ReplacePart[iactive, INPUT_ -> 1] & /@ 
#       Flatten[Position[iactive, 0]];
#     results = (Bestsmap[INPUT_, i, ddbDtg]) & /@ newIactives;
#     mses = results[[All, 1]];
#     mseBest = Min[mses];
#     q = (1 - mseBest/msePrev);
#     hBest = Position[mses, mseBest][[1, 1]];
#     If[q > qc, iactive = newIactives[[hBest]]; msePrev = mseBest; 
#      ciStar = results[[hBest, 2]]; thetaBest = results[[hBest, 3]]];
#     links++
#     ];
#    {ciStar, thetaBest},
#    {baggingIteration}, DistributedContexts -> Automatic, 
#    Method -> "FinestGrained"];
#  {Median[baggingResult[[All, 2]]],
#   Drop[Median[baggingResult[[All, 1]]], -1]}
#  )