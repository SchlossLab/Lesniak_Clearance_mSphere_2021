#!/bin/bash

# Author: Begum Topcuoglu
# Date: 2018-02-13
#
#######################################################################################
# This script will:
#   1. Take single .csv files that have thL2_Logistic_Regressionresult for one datasplit
#   2. Combine them together to have the results for 100 datasplits in one .csv file
#   3. We don't keep the header of each file when combined but only once.

# In the end, the combined_best file must be 101 lines. 1st line is the header and the 100 lines have the data of 100 files.
#             the combined_all file must have 100*(hyper-parameter number)+1 lines.
########################################################################################

SEARCH_DIR=data/temp
FINAL_DIR=data/process

# Keep the first line of File1 and remove the first line of all the others and combine
head -1 $SEARCH_DIR/all_hp_results_L2_Logistic_Regression_1.csv  > $SEARCH_DIR/combined_all_hp_results_L2_Logistic_Regression.csv; tail -n +2 -q $SEARCH_DIR/all_hp_results_L2_Logistic_Regression_*.csv >> $SEARCH_DIR/combined_all_hp_results_L2_Logistic_Regression.csv
head -1 $SEARCH_DIR/best_hp_results_L2_Logistic_Regression_1.csv  > $SEARCH_DIR/combined_best_hp_results_L2_Logistic_Regression.csv; tail -n +2 -q $SEARCH_DIR/best_hp_results_L2_Logistic_Regression_*.csv >> $SEARCH_DIR/combined_best_hp_results_L2_Logistic_Regression.csv
head -1 $SEARCH_DIR/all_imp_features_non_cor_results_L2_Logistic_Regression_1.csv > $SEARCH_DIR/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv; tail -n +2 -q $SEARCH_DIR/all_imp_features_non_cor_results_L2_Logistic_Regression_*.csv >> $SEARCH_DIR/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv
head -1 $SEARCH_DIR/all_imp_features_cor_results_L2_Logistic_Regression_1.csv > $SEARCH_DIR/combined_all_imp_features_cor_results_L2_Logistic_Regression.csv; tail -n +2 -q $SEARCH_DIR/all_imp_features_cor_results_L2_Logistic_Regression_*.csv >> $SEARCH_DIR/combined_all_imp_features_cor_results_L2_Logistic_Regression.csv
head -1 $SEARCH_DIR/walltime_L2_Logistic_Regression_1.csv  > $SEARCH_DIR/walltime_L2_Logistic_Regression.csv; tail -n +2 -q $SEARCH_DIR/walltime_L2_Logistic_Regression_*.csv >> $SEARCH_DIR/walltime_L2_Logistic_Regression.csv
head -1 $SEARCH_DIR/traintime_L2_Logistic_Regression_1.csv  > $SEARCH_DIR/traintime_L2_Logistic_Regression.csv; tail -n +2 -q $SEARCH_DIR/traintime_L2_Logistic_Regression_*.csv >> $SEARCH_DIR/traintime_L2_Logistic_Regression.csv

mv $SEARCH_DIR/traintime_L2_Logistic_Regression.csv $FINAL_DIR/traintime_L2_Logistic_Regression.csv
mv $SEARCH_DIR/combined_all_hp_results_L2_Logistic_Regression.csv $FINAL_DIR/combined_all_hp_results_L2_Logistic_Regression.csv
mv $SEARCH_DIR/combined_best_hp_results_L2_Logistic_Regression.csv $FINAL_DIR/combined_best_hp_results_L2_Logistic_Regression.csv
mv $SEARCH_DIR/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv $FINAL_DIR/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv
mv $SEARCH_DIR/combined_all_imp_features_cor_results_L2_Logistic_Regression.csv $FINAL_DIR/combined_all_imp_features_cor_results_L2_Logistic_Regression.csv

