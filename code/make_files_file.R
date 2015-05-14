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
# * data/mothur/abx_time.files
#
################################################################################

metadata <- read.table(file="data/raw/abx_cdiff_metadata.tsv", header=T)
samples <- rownames(metadata)

#get the list of fastq file names
fastqs <- list.files("../../cdiff_fastqs/", pattern="*.fastq")
have_data <- file.info(paste("../../cdiff_fastqs/", fastqs, sep=""))$size != 0
good_fastqs <- fastqs[have_data]

#need a function that will build the files file line for each stub
find_files <- function(stub, files=good_fastqs){
    pattern <- paste0(stub, ".*R1.*fastq")  #use the stub to create a search pattern
    r1_file <- files[grep(pattern, files)]  #find all of the possible R1 files
    r2_file <- gsub("R1", "R2", r1_file)    #find all of the possible R2 files
    label <- rep(stub, length(r1_file))     #make a vector of the sample name

    #this all is necessary because some samples were sequenced multiple times
    #the output of the above will be a series of vectors that here we will bind
    #as columns, and paste across the columns for each row
    apply(cbind(label, r1_file, r2_file), 1, paste, collapse="\t")
}

files_lines <- lapply(samples, find_files)
write(unlist(files_lines), "data/mothur/abx_time.files")


r1_mock <- good_fastqs[grep("mock[^17].*R1.*fastq", good_fastqs)]
r2_mock <- good_fastqs[grep("mock[^17].*R2.*fastq", good_fastqs)]
stub_mock <- gsub("(.*)_S.*", "\\1", r1_mock)

lines_mock <- apply(cbind(stub_mock, r1_mock, r2_mock), 1, paste, collapse="\t")
write(lines_mock, "data/mothur/abx_time.files", append=T)
