################################################################################
#
#	Part 1: Get the reference files
#
#	Here we give instructions on how to get the necessary reference files that
#	are used throughout the rest of the analysis. These are used to calculate
#	error rates, generate alignments, and classify sequences. We will pull down
#	the mock community reference (HMP_MOCK.fasta), the silva reference alignment
#	(silva.bacteria.align), and the RDP training set data (trainset9_032012).
#	Finally, we use the HMP_MOCK.align to get the alignment coordinates for the
#	V4 data. These data will be stored in the data/references/ folder.
#
#	The targets in this part are all used as dependencies in other rules
#
################################################################################

#Location of files
REFS = data/references/
RAW = data/raw/
MOTHUR = data/mothur/
PROC = data/process/

# create all data and reference files
#get v4 region of the silva reference alignment, rdp training set data and hmp mock
$(REFS)silva.v4.align $(REFS)trainset16_022016.v4.fasta $(REFS)trainset16_022016.v4.tax $(REFS)HMP_MOCK.v4.align : code/get_references.batch
	bash code/get_references.batch

references : $(REFS)silva.v4.align $(REFS)trainset16_022016.v4.fasta $(REFS)trainset16_022016.v4.tax $(REFS)HMP_MOCK.v4.align

################################################################################
#
#	Part 2: Get fastq files
#
################################################################################

# Get metadata from Schubert mBio 2015 doi: 10.1128/mBio.00974-15
data/raw/abxD01_IDS.xlsx : 
	wget -N -P data/raw https://github.com/SchlossLab/Schubert_AbxD01_mBio_2015/raw/master/data/raw/abxD01_IDS.xlsx

# Clean metadata for this project
data/process/abx_cdiff_metadata_clean.txt : data/raw/abxD01_IDS.xlsx \
											code/clean_metadata.R
	Rscript code/clean_metadata.R

# get the fastq files
#### UPDATE USING SRA to download fastqs into data/mothur ####
#$(RAW)get_data : code/get_fastqs.sh $(MOTHUR)abx_clearance.files
#	bash code/get_fastqs.sh $(MOTHUR)abx_clearance.files;\
#	touch $(RAW)get_data

# build the files file. 
$(MOTHUR)abx_clearance.files : code/make_files_file.R $(RAW)abx_cdiff_metadata.tsv
	Rscript code/make_files_file.R


################################################################################
#
#	Part 3: Run data through mothur 
#		Time to run mothur pipline on HPC: 17h:21min
#
################################################################################

SUB = 1200
# here we go from the raw fastq files and the files file to generate a fasta,
# taxonomy, and count_table file that has had the chimeras removed as well as
# any non bacterial sequences
# then we get the sequencing error as seen in the mock community samples
# then we go from the good sequences and generate a shared file and a cons.taxonomy 
# file based on OTU data. Lastly, we rarefy the number of reads to $(SUB) sequences 
# per library for the alpha and beta diversity analyses and modeling and 
# for the shared file
$(MOTHUR)complete.sample.final.shared $(MOTHUR)sample.final.shared $(PROC)mock.sample.final.shared $(MOTHUR)final.taxonomy $(MOTHUR)sample.error.count : code/get_good_seqs.batch\
										code/get_good_seqs.batch\
										code/get_error.batch\
										code/get_shared_otus.batch\
										$(REFS)silva.v4.align\
										$(REFS)trainset16_022016.v4.fasta\
										$(REFS)trainset16_022016.v4.tax\
										$(REFS)HMP_MOCK.v4.fasta
	mothur code/get_good_seqs.batch
	mothur code/get_error.batch
	mothur code/get_shared_otus.batch
	mv data/mothur/abx_clearance.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared data/mothur/complete.sample.final.shared
	mothur "#remove.groups(shared=data/mothur/complete.sample.final.shared, groups=mock2-mock3-mock4-mock5-mock6-mock8-mock9)"
	mv data/mothur/complete.sample.final.0.03.pick.shared data/mothur/sample.final.shared
	mothur "#set.dir(input=data/mothur, output=data/mothur);
		summary.single(shared=sample.final.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=$(SUB));
		dist.shared(shared=sample.final.shared, calc=thetayc, subsample=$(SUB));
		sub.sample(shared=sample.final.shared, size=$(SUB));
		get.groups(shared=complete.sample.final.shared, groups=mock2-mock3-mock4-mock5-mock6-mock8-mock9)"
	mv data/mothur/complete.sample.final.0.03.pick.shared data/process/mock.sample.final.shared
	mv data/mothur/abx_clearance.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy data/mothur/final.taxonomy
	mv data/mothur/abx_clearance.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.error.count data/mothur/sample.error.count


################################################################################
#
#	Part 4: Write the paper
#
################################################################################


################################################################################
# Process raw data for analysis
################################################################################

# Convert taxonomy file into a dataframe with OTU labels
data/process/abx_cdiff_taxonomy_clean.tsv : data/mothur/final.taxonomy\
											code/convert_OTU_labels.R
	Rscript code/convert_OTU_labels.R


################################################################################
#
# Run L2 Logistic Regression
#		Time to run L2 pipline on HPC: 1min
#
################################################################################

SEARCH_DIR=data/temp/l2_otu
FINAL_DIR=data/process/l2_otu
# Create dataframe of subset samples classification column (cleared), and features (OTUs)
# and create correlation matrix of features
data/process/otu_input_data.csv data/process/otu_sample_names.txt data/process/sig_flat_corr_matrix_otu.csv : code/R/setup_model_data.R\
																	code/R/compute_correlation_matrix.R\
																	data/process/abx_cdiff_metadata_clean.txt\
																	data/process/abx_cdiff_taxonomy_clean.tsv\
																	data/mothur/sample.final.0.03.subsample.shared\

	Rscript code/R/setup_model_data.R l2_otu

# Run pipeline array
$SEARCH_DIR/walltime_L2_Logistic_Regression_1.csv : code/R/main.R\
						code/R/run_model.R\
						code/R/model_pipeline.R\
						code/R/tuning_grid.R\
						code/R/permutation_importance.R\
						code/run_logit.sbat\
						code/R/auprc.R\
						code/R/functions.R
	for seed in {1..100}
	do
		Rscript code/R/main.R --seed $seed --model L2_Logistic_Regression --level l2_otu --data  data/process/l2_otu_input_data.csv --hyperparams data/default_hyperparameters.csv --outcome clearance --permutation
	done
	# or run SBATCH code/run_logit.sbat on the Great Lakes cluster

# concatenate results and determine feature importance
$FINAL_DIR/combined_all_imp_features_cor_results_L2_Logistic_Regression.csv : code/bash/process_l2_output.sh\
																				code/R/get_feature_rankings.R
	bash code/bash/process_l2_output.sh


################################################################################
# Create figures
################################################################################

# Figure 1
results/figures/figure_1.jpg : code/build_fig1.R\
								data/process/abx_cdiff_metadata_clean.txt\
								data/mothur/sample.final.0.03.subsample.shared\
								data/process/abx_cdiff_taxonomy_clean.tsv\
								code/sum_otu_by_taxa.R\
								data/mothur/sample.final.groups.ave-std.summary\
								data/mothur/sample.final.thetayc.0.03.lt.ave.dist\
								code/read.dist.R
	Rscript code/build_fig1.R

# Figure 2
results/figures/figure_2.jpg : code/build_fig2.R\
								data/process/abx_cdiff_metadata_clean.txt\
								data/mothur/sample.final.0.03.subsample.shared\
								data/process/abx_cdiff_taxonomy_clean.tsv\
								code/sum_otu_by_taxa.R
	Rscript code/build_fig2.R

# Figure 3
results/figures/figure_3.jpg results/figures/figure_S1.jpg results/figures/figure_S2.jpg : code/build_fig3.R\
								data/process/abx_cdiff_metadata_clean.txt\
								data/mothur/sample.final.0.03.subsample.shared\
								data/process/abx_cdiff_taxonomy_clean.tsv\
								code/sum_otu_by_taxa.R\
								data/mothur/sample.final.groups.ave-std.summary\
								data/mothur/sample.final.thetayc.0.03.lt.ave.dist\
								code/read.dist.R
	Rscript code/build_fig3.R

# Figure S3
results/figures/figure_S3.jpg : code/build_figS3.R\
								data/process/abx_cdiff_metadata_clean.txt\
								data/mothur/sample.final.0.03.subsample.shared\
								data/process/abx_cdiff_taxonomy_clean.tsv\
								code/sum_otu_by_taxa.R\
								data/mothur/sample.final.groups.ave-std.summary\
								data/mothur/sample.final.thetayc.0.03.lt.ave.dist\
								code/read.dist.R
	Rscript code/build_figS3.R

# Figure 4
results/figures/figure_4.jpg results/figures/figure_S4.jpg : code/build_fig4.R\
								code/R/functions.R\
								data/process/abx_cdiff_metadata_clean.txt\
								data/process/abx_cdiff_taxonomy_clean.tsv\
								data/process/otu/combined_best_hp_results_L2_Logistic_Regression.csv\
								data/process/otu/combined_all_sample_results_L2_Logistic_Regression.csv\
								data/process/otu/combined_L2_Logistic_Regression_feature_ranking.tsv\
								data/process/otu/combined_best_hp_results_L2_Logistic_Regression.csv\
								data/process/otu/L2_Logistic_Regression_non_cor_importance.tsv\
								data/process/otu/combined_all_imp_features_non_cor_results_L2_Logistic_Regression.csv\
								data/mothur/sample.final.0.03.subsample.shared
	Rscript code/build_fig4.R otu

# Figure 5
results/figures/figure_5.jpg : code/build_fig5.R\
								data/process/abx_cdiff_metadata_clean.txt\
								data/mothur/sample,final.shared\
								data/mothur/sample.final.0.03.subsample.shared\
								data/process/abx_cdiff_taxonomy_clean.tsv
	Rscript code/build_fig5.R

################################################################################
# commands for processing the dependencies to create the targets

target : dependencies
	commands

write.paper : $(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.cons.taxonomy\
		$(BASIC_STEM).pick.pick.pick.an.unique_list.groups.ave-std.summary\
		$(BASIC_STEM).pick.v4.wang.pick.pick.tx.5.cons.taxonomy\
		$(BASIC_STEM).pick.v4.wang.pick.pick.tx.5.subsample.shared\
		$(BASIC_STEM).pick.pick.pick.error.summary\
		$(MOTHUR)abxD1.counts
