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
save_dir <- paste0('data/process/ccm/interactions')
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
true_interactions <- read.table(paste0('data/process/validation/validation_full_interaction_matrix_', 
  treatment_subset, '.txt')) %>% 
  mutate(affected_otu = otu_list) %>% 
  gather(affector_otu, actual_strength, -affected_otu) %>% 
  mutate(affector_otu =  gsub('V', 'OTU_', affector_otu))
#iters <- 5;i <- 'OTU_3'; j <- c('OTU_6', 'OTU_2', 'OTU_4'); univariate <- F
run_smap <- function(i, j, univariate = F){
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
    make_block(select(composite_ts, i),
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
    if(univariate == T){
      uni_theta <- block_lnlp(select(univariate_ts, contains(i)),
        method = 's-map',
        num_neighbors = 0,
        lib = lib_segments, pred = pred_segments,
        theta = c(0, 1e-04, 3e-04, 0.001,
          0.003, 0.01, 0.03, 0.1,
                  0.3, 0.5, 0.75, 1, 1.5,
                  2, 3, 4, 6, 8),
        target_column = i,
        silent = T)
      test_theta[[iter]] <- data.frame(uni_theta, 
        model_type = 'univariate', stringsAsFactors = F)
      } else {
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
      test_theta[[iter]] <- data.frame(theta_run, 
        model_type = 'multivariate', stringsAsFactors = F)
      }
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
    if(univariate == T){
      smap_uni <-  block_lnlp(select(univariate_ts, contains(i)),
        lib = segments, pred = segments,
        method = "s-map",
        num_neighbors = 0, 
        theta = pull(filter(best_theta, model_type == 'univariate'), theta),
        target_column = i,
        silent = T,
        save_smap_coefficients = T) # save S-map coefficients    
      interaction_smap[[iter]] <- data.frame(smap_uni$model_output[[1]], embed = 'uni', run = iter, 
        mae = smap_uni$mae, stringsAsFactors = F)
    } else {
      feature_names <- colnames(select(rnd_composite_ts, one_of(i, j)))
      smap_multi <-  block_lnlp(select(rnd_composite_ts, one_of(i, j)),
        lib = segments, pred = segments,
        method = "s-map",
        num_neighbors = 0, 
        theta = pull(filter(best_theta, model_type == 'multivariate'), theta),
        target_column = i, columns = c(i, j),
        silent = T,
        save_smap_coefficients = T) # save S-map coefficients
      smap_coef_df <- smap_multi$smap_coefficients[[1]]
      colnames(smap_coef_df) <- c(paste0('d', i, '_d', feature_names), 'intercept')
      interaction_smap[[iter]] <- data.frame(bind_cols(smap_multi$model_output[[1]], smap_coef_df,
          right_join(composite_ts, pred_order, by = 'unique_id')), 
        embed = 'multi', run = iter, mae = smap_multi$mae, stringsAsFactors = F)
    }
  }
  output <- list(model_output = do.call('rbind', interaction_smap))
  if(univariate == T){
    output <- data.frame(affected_otu = i, 
      median_mae = median(output[['model_output']]$mae), stringsAsFactors = F)
  } else {
    output <- full_join(
        output[['model_output']] %>% 
            group_by(embed) %>%
            summarise(median_mae = median(mae)),
        output[['model_output']] %>% 
          gather(interaction, strength, contains('dOTU')) %>%
          group_by(interaction, embed) %>% 
          summarise(median_interaction = median(strength, na.rm = T),
            ixn_lower_qrtl = quantile(strength, na.rm = T)['25%'],
            ixn_upper_qrtl = quantile(strength, na.rm = T)['75%']),
          by = 'embed') %>% 
          mutate(otus_tmp = gsub('d', '', gsub('_d', 'SEP', interaction))) %>% 
      separate(otus_tmp, c('affected_otu', 'affector_otu'), sep = 'SEP')
  }
  return(output)
}

Qc <- 0.01

# check univariate performance and set as initial MSE
univariate_smap_output <- map_dfr(otu_list, function(i) run_smap(i, i, univariate = T))
MSE_prev <- univariate_smap_output

# for each otu
#iters <- 500
# run smap stepwise and incorporate the best OTU at each step
  # run all otus by all other OTUs, keeping the best performing 
  # and then checking for improvement by remaining OTUs 
  # and repeat until additonal OTU doesnt provide improvement
smap_output <- map_dfr(otu_list, function(current_otu){
  Q <- 1
  remaining_otus <- otu_list[!otu_list %in% current_otu]
  active <- c()
  MSE_current <- MSE_prev[MSE_prev$affected_otu == current_otu, 'median_mae']
  while(Q > Qc & length(active$affector_otu) < length(remaining_otus)){
    # create list of predictive OTUs with one of each OTU not most predictive yet
    stepwise_otu_list <- lapply(as.list(otu_list[!otu_list %in% c(current_otu, active$affector_otu)]),
      function(x) append(x, active$affector_otu))
    # run smap with set of predictive OTUs with one new OTU (for all unactive OTUs)
    smap_out <- map_dfr(stepwise_otu_list, function(j) run_smap(current_otu, j))
    # identify which of the new OTUs results in the best performance (lowest MAE)
    MSE_best <- smap_out %>% 
    	filter(!affector_otu %in% c(current_otu, active$affector_otu)) %>%
     	filter(median_mae == min(median_mae)) %>% 
      select(affector_otu, median_mae)
    # determine the relative increase in MAE
    Q <- 1 - (MSE_best$median_mae / MSE_current)
    # if increase in MAE is greater than threshold (Qc) then add OTU to active and update MSE
    if(Q > Qc){
      active <- rbind(active, data.frame(MSE_best, Q = Q))
      MSE_current <- MSE_best$median_mae
      best_interaction <- smap_out %>% 
        filter(median_mae == min(median_mae)) %>% 
        select(affector_otu, affected_otu, median_interaction, ixn_upper_qrtl, ixn_lower_qrtl)
      }
  }
  if(length(active) > 0){
      return(full_join(rbind(data.frame(affector_otu = current_otu, affected_otu = current_otu, 
            median_mae = MSE_prev[MSE_prev$affected_otu == current_otu, 'median_mae'], Q = 0, stringsAsFactors = F),
          data.frame(affected_otu = current_otu, active, stringsAsFactors = F)),
        best_interaction, by = c('affector_otu', 'affected_otu')))
    } else {
      return(full_join(rbind(data.frame(affector_otu = current_otu, affected_otu = current_otu, 
            median_mae = MSE_prev[MSE_prev$affected_otu == current_otu, 'median_mae'], Q = 0, stringsAsFactors = F),
          data.frame(affected_otu = current_otu, affector_otu = 'NA', median_mae = NA, Q = NA, stringsAsFactors = F)),
        best_interaction, by = c('affector_otu', 'affected_otu')))

    }
})

smap_output <- smap_output  %>% 
  group_by(affected_otu) %>% 
  mutate(order = (1:length(affector_otu)) - 1) %>% 
  ungroup

detected_ixn_df <- full_join(true_interactions, smap_output, by = c('affector_otu', 'affected_otu')) %>% 
  mutate(actual_ixn_TF = ifelse(abs(actual_strength) > 0, T, F),
    predicted_ixn_TF = ifelse(is.na(median_mae), F, T))

write.table(detected_ixn_df, paste0(save_dir, '/interactions_', treatment_subset, '.txt'), 
  quote = F, row.names = F)


Actual_pos <- nrow(filter(detected_ixn_df, actual_ixn_TF == T)) 
TP <- nrow(filter(detected_ixn_df, actual_ixn_TF == T & predicted_ixn_TF == T)) 
sensitivity <- nrow(filter(detected_ixn_df, actual_ixn_TF == T & predicted_ixn_TF == T)) / 
  nrow(filter(detected_ixn_df, actual_ixn_TF == T))
specificity <- nrow(filter(detected_ixn_df, actual_ixn_TF == F & predicted_ixn_TF == F)) / 
  nrow(filter(detected_ixn_df, actual_ixn_TF == F))

data.frame(treatment = treatment_subset, sensitivity = sensitivity, specificity = specificity,
  detected_ixn_df)



##   # comparison to univariate is limited to late time points 
##   # so  limited interpretation
##   #interaction_smap %>% 
##   # ggplot(aes(x = obs, y = pred, color = embed)) + 
##   #   geom_point(alpha = .2) + 
##   #   theme_bw() + 
##   #     coord_cartesian(ylim=c(-2, 2), xlim=c(-2, 2))
##   #interaction_smap %>% 
##   # group_by(embed) %>% 
##   # summarise(mean = mean(obs, na.rm = T), median = median(obs, na.rm = T),
##   #   pval = cor.test(.$obs, .$pred, method = 'spearman')$p.val,
##   #   cor = cor.test(.$obs, .$pred, method = 'spearman')$estimate)
## 
##   ## Time series of fluctuating interaction strength
##   order <- expand.grid(unique_id = sample_list, day = 0:10)  %>% 
##     arrange(unique_id, day) %>% 
##     mutate(epochs = 1:NROW(.)) #%>% 
##     #right_join(expand.grid(epochs = 1:NROW(.), run = 1:500)) %>% 
##     #arrange(unique_id, run, day)
##   # Plot all partial derivatives
##   interaction_temporal_plot <- interaction_smap %>% 
##     filter(embed == 'multi') %>% 
##     right_join(order, by = c('day', 'unique_id')) %>% 
##     gather(interaction, strength, one_of(paste0('d', i, '_d', j))) %>% 
##     mutate(interaction = gsub('tu0*', 'TU', interaction)) %>% 
##     ggplot(aes(x = epochs, y = strength)) + 
##       stat_summary(fun.data = 'median_hilow', geom = 'ribbon', 
##         alpha = 0.2, fun.args =(conf.int = 0.5), aes(fill = interaction)) + 
##       stat_summary(aes(color = interaction), fun.y = median, geom = 'line') + 
##       theme_bw() +
##       theme(legend.position="top", legend.title=element_blank()) +
##       labs(title = paste0('Strength of effect of OTUs on ', gsub('tu0*', 'TU', i)),
##         subtitle = 'Median with IQR shown')
## 
##   interaction_dist_plot <- interaction_smap %>% 
##     filter(embed == 'multi') %>% 
##     right_join(order, by = c('day', 'unique_id')) %>% 
##     gather(interaction, strength, one_of(paste0('d', i, '_d', j))) %>% 
##     mutate(interaction = gsub('^.*_d', '', interaction)) %>% 
##     inner_join(filter(true_interactions, affected_otu == i), ) %>% 
##     ggplot(aes(x = interaction, y = strength, fill = interaction)) + 
##       geom_hline(yintercept = 0, linetype = 'dashed') + 
##       geom_violin(alpha = 0.5) + 
##       geom_boxplot(aes(y = actual_strength), color = 'red') + 
##       theme_bw() +
##       theme(legend.position='none') +
##       labs(title = paste0('Strength of effect of OTUs on ', gsub('tu0*', 'TU', i)),
##         subtitle = 'Median with IQR shown')
##   
## # dynamics_plot <- shared_file %>% 
## #     select(cage, mouse, day, one_of(i, j)) %>% 
## #     mutate(unique_id = paste(cage, mouse, sep = '_')) %>% 
## #     right_join(order) %>% 
## #     #mutate(time = 1:nrow(.)) %>% 
## #     gather(bacteria, counts, one_of(i, j)) %>% 
## #       ggplot(aes(x = epochs, y = counts, color = as.factor(cage), group = interaction(cage, mouse))) + 
## #         geom_line(alpha = 0.4) + 
## #         geom_point() + 
## #         facet_grid(bacteria~., scales = 'free_y') +
## #         theme_bw() + 
## #         labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
## #           subtitle = 'Temporal Dynamics - Colored by mouse', 
## #           title = paste(i, 'is driven by', paste(j, collapse = ', '))) + 
## #         theme_bw(base_size = 8) + 
## #         theme(legend.position = 'none', panel.grid.minor = element_blank())
## 
##   ggsave(filename = paste0(save_dir, '/interactions_w_', i, '.jpg'), 
##     plot = plot_grid(interaction_temporal_plot, interaction_dist_plot, nrow = 2, align = 'v'),
##     width = 14, height = 14, device = 'jpeg')
##   
## 
## #}
## }
## # output file with taxa, interaction direction, interaction strength
## 


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