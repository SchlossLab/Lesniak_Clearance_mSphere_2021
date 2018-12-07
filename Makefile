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

#get the silva reference alignment
$(REFS)silva.bacteria.align :
	wget -N -P $(REFS) http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz; \
	tar xvzf $(REFS)Silva.nr_v119.tgz -C $(REFS);
	mothur "#get.lineage(fasta=$(REFS)silva.nr_v119.align, taxonomy=$(REFS)silva.nr_v119.tax, taxon=Bacteria)";
	mv $(REFS)silva.nr_v119.pick.align $(REFS)silva.bacteria.align; \
	rm $(REFS)README.html; \
	rm $(REFS)README.Rmd; \
	rm $(REFS)silva.nr_v119.*

#get the v4 region of the alignment
$(REFS)silva.v4.align : $(REFS)silva.bacteria.align
	mothur "#pcr.seqs(fasta=$(REFS)silva.bacteria.align, start=11894, end=25319, keepdots=F, processors=8);\
			unique.seqs(fasta=current);"; \
	mv $(REFS)silva.bacteria.pcr.unique.align $(REFS)silva.v4.align; \
	rm $(REFS)silva.bacteria.pcr.*

#get the rdp training set data
$(REFS)trainset10_082014.pds.tax $(REFS)trainset10_082014.pds.fasta :
	wget -N -P $(REFS) http://www.mothur.org/w/images/2/24/Trainset10_082014.pds.tgz; \
	tar xvzf $(REFS)Trainset10_082014.pds.tgz -C $(REFS);\
	mv $(REFS)trainset10_082014.pds/trainset10_082014.* $(REFS);\
	rm -rf $(REFS)trainset10_082014.pds

#get the V4 region of the RDP training set
$(REFS)trainset10_082014.v4.tax $(REFS)trainset10_082014.v4.fasta : \
						$(REFS)trainset10_082014.pds.tax \
						$(REFS)trainset10_082014.pds.fasta \
						$(REFS)silva.v4.align
	mothur "#align.seqs(fasta=$(REFS)trainset10_082014.pds.fasta, reference=$(REFS)silva.v4.align, processors=8);\
		screen.seqs(fasta=current, taxonomy=$(REFS)trainset10_082014.pds.tax, start=1, end=13425);\
		degap.seqs(fasta=current)"; \
	mv $(REFS)trainset10_082014.pds.good.ng.fasta $(REFS)trainset10_082014.v4.fasta; \
	mv $(REFS)trainset10_082014.pds.good.tax $(REFS)trainset10_082014.v4.tax;\
	rm $(REFS)trainset10_082014.pds.align*;\
	rm $(REFS)trainset10_082014.pds.bad.accnos;\
	rm $(REFS)trainset10_082014.pds.flip.accnos;

$(REFS)HMP_MOCK.fasta :
	wget --no-check-certificate -N -P $(REFS) https://raw.githubusercontent.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/master/data/references/HMP_MOCK.fasta

#align the mock community reference sequeces
$(REFS)HMP_MOCK.v4.fasta : $(REFS)HMP_MOCK.fasta $(REFS)silva.v4.align
	mothur "#align.seqs(fasta=$(REFS)HMP_MOCK.fasta, reference=$(REFS)silva.v4.align);\
			degap.seqs()";\
	mv $(REFS)HMP_MOCK.ng.fasta $(REFS)HMP_MOCK.v4.fasta;\
	rm $(REFS)HMP_MOCK.align;\
	rm $(REFS)HMP_MOCK.align.report;\
	rm $(REFS)HMP_MOCK.flip.accnos

references : $(REFS)HMP_MOCK.v4.fasta $(REFS)trainset10_082014.v4.tax $(REFS)trainset10_082014.v4.fasta $(REFS)silva.v4.align

################################################################################
#
#	Part 2: Get fastq files
#
################################################################################


# build the files file. probably should replace this chunk eventually
# with pulling data off of the SRA
#$(MOTHUR)abx_time.files : code/make_files_file.R $(RAW)abx_cdiff_metadata.tsv
#	R -e "source('code/make_files_file.R')"


# need to get the fastq files. probably should replace this chunk eventually
# with pulling data off of the SRA
#$(RAW)get_data : code/get_fastqs.sh $(MOTHUR)abx_time.files
#	bash code/get_fastqs.sh $(MOTHUR)abx_time.files;\
#	touch $(RAW)get_data

# need to get the CFU on the day after antibiotic treatment along with the
# part of the experiment that each sample belongs to
#
#$(MOTHUR)abxD1.counts : code/make_counts_file.R $(MOTHUR)abx_time.files\
#							$(MOTHUR)abx_cdiff_metadata.tsv
#	R -e "source('code/make_counts_file.R')"



################################################################################
#
#	Part 3: Run data through mothur
#
################################################################################

BASIC_STEM = $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster

# here we go from the raw fastq files and the files file to generate a fasta,
# taxonomy, and count_table file that has had the chimeras removed as well as
# any non bacterial sequences
$(BASIC_STEM).uchime.pick.pick.count_table $(BASIC_STEM).pick.pick.fasta $(BASIC_STEM).pick.v4.wang.pick.taxonomy : code/get_good_seqs.batch\
										$(RAW)get_data\
										$(REFS)silva.v4.align\
										$(REFS)trainset10_082014.v4.fasta\
										$(REFS)trainset10_082014.v4.tax
	mothur code/get_good_seqs.batch;\
	rm $(MOTHUR)*.map



# here we go from the good sequences and generate a shared file and a
# cons.taxonomy file based on OTU data
$(BASIC_STEM).pick.pick.pick.an.unique_list.shared $(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.cons.taxonomy : code/get_shared_otus.batch\
										$(BASIC_STEM).uchime.pick.pick.count_table\
										$(BASIC_STEM).pick.pick.fasta\
										$(BASIC_STEM).pick.v4.wang.pick.taxonomy
	mothur code/get_shared_otus.batch;\
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.pick.count_table;\
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta;\
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.taxonomy



# here we go from the good sequences and generate a shared file and a
# cons.taxonomy file based on phylum-level data
$(BASIC_STEM).genus.cons.taxonomy $(BASIC_STEM).phylum.cons.taxonomy $(BASIC_STEM).genus.shared $(BASIC_STEM).phyla.shared : code/get_shared_phylotypes.batch\
										$(BASIC_STEM).uchime.pick.pick.count_table\
										$(BASIC_STEM).pick.pick.fasta\
										$(BASIC_STEM).pick.v4.wang.pick.taxonomy
	mothur code/get_shared_phylotypes.batch;
	grep "^1" $(BASIC_STEM).pick.v4.wang.pick.pick.tx.shared > $(BASIC_STEM).genus.shared
	grep "^5" $(BASIC_STEM).pick.v4.wang.pick.pick.tx.shared > $(BASIC_STEM).phylum.shared
	mv $(BASIC_STEM).pick.v4.wang.pick.pick.tx.1.cons.taxonomy $(BASIC_STEM).genus.cons.taxonomy
	mv $(BASIC_STEM).pick.v4.wang.pick.pick.tx.5.cons.taxonomy $(BASIC_STEM).phylum.cons.taxonomy
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.*.cons.tax.summary
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.cons.tax.summary
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.pick.count_table;\
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta;\
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.taxonomy;\
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.equalized.tree.sum
	rm $(MOTHUR)abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.equalized.taxonomy
	rm $(MOTHUR)*.tx.*rabund;


# now we want to get the sequencing error as seen in the mock community samples
$(BASIC_STEM).pick.pick.pick.error.summary : code/get_error.batch\
										$(BASIC_STEM).uchime.pick.pick.count_table\
										$(BASIC_STEM).pick.pick.fasta\
										$(REFS)HMP_MOCK.v4.fasta
	mothur code/get_error.batch


SUB = 2000
# rarefy the number of reads to $(SUB) sequences per library for the alpha and beta diversity analyses and modeling
$(BASIC_STEM).pick.pick.pick.an.unique_list.0.03.subsample.shared $(BASIC_STEM).pick.pick.pick.an.unique_list.groups.ave-std.summary $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.dist : $(BASIC_STEM).pick.pick.pick.an.unique_list.shared
	mothur "#dist.shared(shared=$^, calc=thetayc, subsample=$(SUB), iters=100); summary.single(shared=$^, subsample=$(SUB), calc=nseqs-sobs-shannon-invsimpson, iters=100); sub.sample(shared=$^, size=$(SUB))";\
	rm $(BASIC_STEM).pick.pick.pick.an.unique_list.groups.summary;\
	rm $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.dist;\
	rm $(BASIC_STEM).pick.pick.pick.an.unique_list.thetayc.0.03.lt.std.dist;\
	rm $(BASIC_STEM).pick.pick.pick.an.unique_list.*.rabund

# rarefy the number of reads to $(SUB) sequences per library for the barcarts
$(BASIC_STEM).genus.1.subsample.shared $(BASIC_STEM).phyla.5.subsample.shared : $(BASIC_STEM).genus.shared $(BASIC_STEM).phylum.shared
	mothur "#sub.sample(shared=$(BASIC_STEM).genus.shared, size=$(SUB))";
	mothur "#sub.sample(shared=$(BASIC_STEM).phylum.shared, size=$(SUB))";


################################################################################
#
#	Part 4: Write the paper
#
################################################################################

# Generate simulated community with known interactions
data/process/validation/* : code/generate_validation_data.R
	for gamma in `seq 0 2`;
	do
		for seed in `seq 1 5`;
		do
		Rscript code/generate_validation_data.R $seed $gamma
		done
	done 

################################################################################
# Setup data for CCM
# Z-score normalize and first-difference data
data/process/abx_cdiff_metadata_clean.txt : code/clean_metadata.R\
											data/raw/abx_cdiff_metadata.tsv 
	Rscript code/clean_metadata.R
data/process/ccm_otu_data.txt : code/setup_ccm_data.R\
								data/process/abx_cdiff_metadata_clean.txt
	Rscript code/setup_ccm_data.R
data/process/validation/ccm_validation_data.txt : data/process/validation/validation_temporal_data_*.txt \
												  data/process/validation/validation_interaction_matrix_*.txt\
												  code/setup_validation_data.R
	Rscript code/setup_validation_data.R

################################################################################
# determine best embedding for each OTU
data/process/ccm/*/simplex_embedding_first_differenced.txt : code/get_best_embedding.R\ 
								data/process/bucci/ccm_bucci_data.txt\
								data/process/ccm_otu_data.txt
	for i in `seq 1 15`;
	do
		Rscript code/get_best_embedding.R $i 'data/process/ccm_otu_data.txt'
	done  
data/process/ccm/validation/simplex_embedding_first_differenced.txt : code/get_best_embedding.R\ 
								data/process/ccm_validation_data.txt
	for i in `seq 1 15`;
	do
		Rscript code/get_best_embedding.R $i 'data/process/ccm_validation_data.txt'
	done  

################################################################################
# check each OTU for nonlinearity
data/process/ccm/*/smap_nonlinearity_first_differenced.txt : code/run_smap.R\ 
								data/process/ccm/*/simplex_embedding_first_differenced.txt
	for i in `seq 1 15`;
	do
		Rscript code/run_smap.R $i
	done  
data/process/ccm/*/smap_nonlinearity_first_differenced.txt : code/run_smap.R\ 
							   data/process/ccm/*/simplex_embedding_first_differenced.txt
	for i in `seq 1 15`;
	do
		Rscript code/run_smap.R $i 'data/process/ccm_validation_data.txt'
	done  

################################################################################
# run CCM on all OTUs that are significantly nonlinear
data/process/ccm/*/ccm_by_otu_*_first_differenced.txt : code/ccm_analysis.R\ 
														data/process/ccm/*/simplex_embedding_first_differenced.txt\
														data/process/ccm/*/smap_nonlinearity_first_differenced.txt
	for i in `seq 1 15`;
	do
		Rscript code/ccm_analysis.R $i 'data/process/ccm_otu_data.txt'
	done  
data/process/ccm/*/ccm_by_otu_*_first_differenced.txt : code/run_smap.R\ 
														data/process/ccm/*/simplex_embedding_first_differenced.txt\
														data/process/ccm/*/smap_nonlinearity_first_differenced.txt
	for i in `seq 1 15`;
	do
		Rscript code/ccm_analysis.R $i 'data/process/ccm_validation_data.txt'
	done  

################################################################################
# run Smap stepwise to detect OTUs most predictive

data/process/ccm/interactions/interactions_*.txt : code/stepwise_smap.R\ 
												   data/process/ccm/ccm_validation_data.txt
	for i in `seq 1 15`;
	do
		Rscript code/stepwise_smap.R $i 'data/process/ccm_validation_data.txt'
	done  

################################################################################
# run smap to determine interactions

data/process/ccm/*/interactions_w_*.txt : code/get_ccm_interactions.R\
										  data/process/ccm_validation_data.txt\
										  data/process/ccm/*/simplex_embedding_first_differenced.txt\
										  data/process/ccm/*/smap_nonlinearity_first_differenced.txt\
										  data/process/ccm/*/ccm_by_otu_*_first_differenced.txt
	Rscript code/get_ccm_interactions.R $i 'data/process/ccm_validation_data.txt'


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
