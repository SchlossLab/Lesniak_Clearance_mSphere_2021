################################################################################
#
# make_files_file.R
#
# This script will build the files file that will be used in the make.contigs
# command from within mothur.
#
# Dependencies...
# * data/mothur/abx_cdiff_metadata.tsv
# * fastq files stored in ../../cdiff_fastqs/
#
# Output...
# * data/mothur/abx_clearance.files
#
################################################################################

library(tidyverse)

# read in meta data to get the group name for each sample
metadata <- read_tsv('data/process/abx_cdiff_metadata_clean.txt',
		col_types = c('cddddcdllcdddddcddcc')) %>% 
# create a sample column with to use to match to the fastq file name
	mutate(sample = paste0(cage, '_', mouse, '_D', day)) %>% 
	select(group, sample, mouse_id)

# get the list of fastq file names
fastqs <- list.files("data/mothur/", pattern="*.fastq$")
have_data <- file.info(paste("data/mothur/", fastqs, sep=""))$size != 0
good_fastqs <- fastqs[have_data]

# create a sample column using the samples cage, mouse ear tag and day to match meta
seqs_by_cage <- tibble(seq_id = fastqs) %>% 
	filter(!grepl('mock', seq_id)) %>% # remove mock
	mutate(seq_R = case_when(grepl('R1', seq_id) ~ 'R1', 
			grepl('R2', seq_id) ~ 'R2',
			T ~ 'NA')) %>% # create column for fastq pairs
	mutate(group = gsub('\\_S\\d.*', '', seq_id)) %>% # remove sequencing filename suffix
	separate(group, c('mouse_id', 'day'), sep = 'D') %>% # separate day
	separate(mouse_id, c('cage', 'mouse'), sep = '-') %>% # separate cage and mouse ids
	mutate(# convert to numeric to remove leading 0s
		sample = paste0(as.numeric(cage), '_', as.numeric(mouse), '_D', as.numeric(day)),
		unique_seq = gsub('\\_R.*', '', seq_id)) %>% # create unique id for samples sequenced more than one
	pivot_wider(names_from = 'seq_R', values_from = 'seq_id') %>% # spread fastqs into R1/R2 columns
	select(sample, R1, R2)
# setup mock fastqs
mock_seqs <- tibble(seq_id = fastqs) %>% 
	filter(grepl('mock[^17]', seq_id)) %>% # remove mock 1 and 7 (as previously done with mothur pipeline originally processing these samples)
	mutate(group = gsub('\\_S\\d.*', '', seq_id),
		seq_R = ifelse(grepl('R1', seq_id), 'R1', 'R2')) %>% 
	pivot_wider(names_from = 'seq_R', values_from = 'seq_id')

# join metadata group ids to fastq file names by created sample id and then join mock
files_df <- metadata %>% 
	inner_join(seqs_by_cage, by = c('sample')) %>% 
	select(group, R1, R2) %>% 
	bind_rows(mock_seqs)

write_tsv(files_df, "data/mothur/abx_clearance.files", col_names = F)
