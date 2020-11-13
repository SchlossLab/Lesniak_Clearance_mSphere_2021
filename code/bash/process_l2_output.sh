#!/bin/bash

# Initial Author: Begum Topcuoglu
# Date: 2018-02-13
# Revised: Nick Lesniak
#
#######################################################################################
# This script will:
#   1. Take single .csv files that have the model result for one datasplit
#   2. Combine them together to have the results for 100 datasplits in one .csv file
#   3. We don't keep the header of each file when combined but only once.

# In the end, the combined_best file must be 101 lines. 1st line is the header and the 100 lines have the data of 100 files.
#             the combined_all file must have 100*(hyper-parameter number)+1 lines.
########################################################################################

for dir in "l2_otu"
do
    SEARCH_DIR=data/temp/$dir
    FINAL_DIR=data/process/$dir
    mkdir $FINAL_DIR
    # Keep the first line of File1 and remove the first line of all the others and combine

    for model in "L2_Logistic_Regression" #"Random_Forest" #"Decision_Tree"
    do
      	head -1 $SEARCH_DIR/all_sample_results_"$model"_1.csv  > $SEARCH_DIR/combined_all_sample_results_"$model".csv; tail -n +2 -q $SEARCH_DIR/all_sample_results_"$model"_*.csv >> $SEARCH_DIR/combined_all_sample_results_"$model".csv
        head -1 $SEARCH_DIR/all_hp_results_"$model"_1.csv  > $SEARCH_DIR/combined_all_hp_results_"$model".csv; tail -n +2 -q $SEARCH_DIR/all_hp_results_"$model"_*.csv >> $SEARCH_DIR/combined_all_hp_results_"$model".csv
        head -1 $SEARCH_DIR/best_hp_results_"$model"_1.csv  > $SEARCH_DIR/combined_best_hp_results_"$model".csv; tail -n +2 -q $SEARCH_DIR/best_hp_results_"$model"_*.csv >> $SEARCH_DIR/combined_best_hp_results_"$model".csv
        head -1 $SEARCH_DIR/all_imp_features_non_cor_results_"$model"_1.csv > $SEARCH_DIR/combined_all_imp_features_non_cor_results_"$model".csv; tail -n +2 -q $SEARCH_DIR/all_imp_features_non_cor_results_"$model"_*.csv >> $SEARCH_DIR/combined_all_imp_features_non_cor_results_"$model".csv
    	head -1 $SEARCH_DIR/all_imp_features_cor_results_"$model"_1.csv > $SEARCH_DIR/combined_all_imp_features_cor_results_"$model".csv; tail -n +2 -q $SEARCH_DIR/all_imp_features_cor_results_"$model"_*.csv >> $SEARCH_DIR/combined_all_imp_features_cor_results_"$model".csv
        head -1 $SEARCH_DIR/walltime_"$model"_1.csv  > $SEARCH_DIR/walltime_"$model".csv; tail -n +2 -q $SEARCH_DIR/walltime_"$model"_*.csv >> $SEARCH_DIR/walltime_"$model".csv
        head -1 $SEARCH_DIR/traintime_"$model"_1.csv  > $SEARCH_DIR/traintime_"$model".csv; tail -n +2 -q $SEARCH_DIR/traintime_"$model"_*.csv >> $SEARCH_DIR/traintime_"$model".csv

        mv $SEARCH_DIR/traintime_"$model".csv $FINAL_DIR/traintime_"$model".csv
        mv $SEARCH_DIR/combined_all_sample_results_"$model".csv $FINAL_DIR/combined_all_sample_results_"$model".csv
        mv $SEARCH_DIR/combined_all_hp_results_"$model".csv $FINAL_DIR/combined_all_hp_results_"$model".csv
        mv $SEARCH_DIR/combined_best_hp_results_"$model".csv $FINAL_DIR/combined_best_hp_results_"$model".csv
    	mv $SEARCH_DIR/combined_all_imp_features_non_cor_results_"$model".csv $FINAL_DIR/combined_all_imp_features_non_cor_results_"$model".csv
        mv $SEARCH_DIR/combined_all_imp_features_cor_results_"$model".csv $FINAL_DIR/combined_all_imp_features_cor_results_"$model".csv
    done
    # get rankings
    Rscript code/R/get_feature_rankings.R $dir

    RANK=feature_ranking

    # 1. Keep the first line of File0 and remove the first line of all the other files (File[0-99]) and
    #       output it to the FINAL_DIR location
    cp $SEARCH_DIR/L2_Logistic_Regression_"$RANK"_1.tsv $FINAL_DIR/combined_L2_Logistic_Regression_"$RANK".tsv


    #   2. Append the other files to the end, but we want to be sure to ignore the 0 file since we don't
    #       want it printed twice
    #        "tail -n +2" makes tail print lines from 2nd line to the end
    #        "-q" tells it to not print the header with the file name
    #        ">>" adds all the tail stuff from every file to the combined file
    tail -n +2 -q $SEARCH_DIR/L2_Logistic_Regression_"$RANK"_{2..100}.tsv >> $FINAL_DIR/combined_L2_Logistic_Regression_"$RANK".tsv
done