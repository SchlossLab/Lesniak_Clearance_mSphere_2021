# Exploring Alyx's time series data, day 1

# to do things interactively on axiom, do the following
# cloned repo to kaitlin on axiom, used link command to copy data over
# navigate to mothur directory

# Run R by typing "R" into the terminal 

# OR- because that's a pain in the ass, scp the important files onto local drive

# scp kjflynn@axiom.ccmb.med.umich.edu:/mnt/EXT/Schloss-data/kaitlin/schuberttime/Schubert_time_XXXX_2015/data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared /Users/Kaitlin/Documents/Schloss_Lab/Flynn_Cdiff_Dynamics/

# Start exploring by recreating some figures from Alyx's 2015 mBio paper (this code is primarily theirs)

# This script will build Figure 1 from their paper, which is a barchart of median relative abundance of each genus found in mice treated with antibiotics as well as untreated control (depending upon what metadata I have)

#read in taxonomy OTU data
taxonomy_file <- read.table(file="abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)

#take taxonomy column full of taxonomy names and assign it to taxonomy
taxonomy <- taxonomy_file$Taxonomy

#assign OTU row names to 'names'
names(taxonomy) <- rownames(taxonomy_file)

#use gsub as regex command to trim taxonomy OTUs to highest level of classification possible
#and remove punctuation

taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub(";unclassified", "", taxonomy)
taxonomy <- gsub("/.*", "", taxonomy)
taxonomy <- gsub(";$", "", taxonomy)
taxonomy <- gsub(".*;", "", taxonomy)

#read in metadata file
metadata <- read.table(file="abx_cdiff_metadata.tsv", header = T)

#DONT DO I THINK... pull out mice that were treated with clinda or control using indexing
#clinda <- metadata[metadata$abx=="clinda" | metadata$abx=="none",]

#read in shared OTU counts file to get relative abundances
shared_file <- read.table(file="abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared", header = T, row.names=2)

#this takes out the 'label' and 'numOTUs' columns from the shared file
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]

#apply sum function to rows (file, 1 for rows, sum for sum). what are brackets? first row?
n_seqs <- apply(shared_file, 1, sum)[1]

#calculate relative abundance percentage for all OTUs in matrix
rel_abund <- 100*shared_file/n_seqs

#find the overlap of samples in rel_abund that come from metadata and organize them by mouse and day
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(metadata))]

#assign overlap samples to row names in rel_abund
rel_abund <- rel_abund[overlap,]

#same idea here if subsampling
#clinda <- clinda[overlap,]

#otherwise
metadata <- metadata[overlap,]

#limit analysis to OTUs that have mean relative abundance over 1% within each antibiotic dose (i feel like i'm getting off track here)
control_rabund <- rel_abund[metadata$abx == "none",]
#get metadata for control mice
control_metadata <- metadata[metadata$abx == "none",]
#find median over columns
control_median <- apply(control_rabund, 2, median)
#returns logical of if OTUs are above 1 percent
control_otus <- control_median > 1.0

#Now we are going to use a function to limit analysis for the drug datasets

otu_hyp_test <- function(drug){

#limit analysis to OTUs that have mean RA over 3%
    drug_rabund <- rel_abund[metadata$abx == drug,]
    drug_metadata <- metadata[metadata$abx == drug,]
    drug_otus <- apply(drug_rabund, 2, median) > 3.0

    combined_otus <- control_otus | drug_otus
    #bind the drug and control into a new vector, use row names, then combined otu are cols?
    combined_rabund <- rbind(drug_rabund, control_rabund)[combined_otus]

    #building this new vector by adding a drugged colum for t/f
    drugged <- c(rep(TRUE, nrow(drug_rabund)), rep(FALSE, nrow(control_rabund)))
    n_otus <- ncol(combined_rabund)

    #make empty p_value column
    p_values <- rep(NA, n_otus)

    #set options to ignore all warning messages ("warn"= -1)
    warn_orig <- options("warn")$warn
    options("warn"= -1)

    #loop over empty p value column and do wilcoxon test on corresponding data column
    for(i in 1:n_otus){
        p_values[i] <- wilcox.test(combined_rabund[,i], g=drugged)$p.value
    }

    #set warning back?
    options("warn" = warn_orig)
    #adjust p values using BH correction and assign to column
    adj_p_values <- p.adjust(p_values, method="BH")
    colnames(combined_rabund)[adj_p_values<0.05]
}

#run function on each individual abx sample, store in variable
amp_sig_otus <- otu_hyp_test("amp")
cef_sig_otus <- otu_hyp_test("cef")
cipro_sig_otus <- otu_hyp_test("cipro")
clinda_sig_otus <- otu_hyp_test("clinda")
metro_sig_otus <- otu_hyp_test("metro")
strep_sig_otus <- otu_hyp_test("strep")
vanc_sig_otus <- otu_hyp_test("vanc")

#now sort, extracting unique elements from the list of otus 
sig_otus <- sort(unique(c(amp_sig_otus, cef_sig_otus, cipro_sig_otus,
                          clinda_sig_otus, metro_sig_otus, strep_sig_otus,
                          vanc_sig_otus)))

#not entirely sure why this exists but it multiplies the length by 1.2
x_max <- length(sig_otus) * 1.2

#also mostly unsure here. says what order each OTU is in rel abundance for control median?. so first listed is 7th abundant
#but only in first row???
o <- order(control_median[sig_otus], decreasing=T)

#organizes relative abundance by significant otus and o order but order only correct for first mouse?
rel_abund_sig <- rel_abund[,sig_otus[o]]

# change labels from OTU000 to taxonomy name 
otu <- gsub("Otu0*", "OTU~", names(taxonomy))
names(otu) <- names(taxonomy)


#put OTU and taxonomy name is one column/variable  
tax_label <- paste0("italic('", taxonomy[sig_otus[o]], "')~(", otu[sig_otus[o]], ")")

# set up function to plot bars of data for drugs
single_drug_bars <- function(drug, drug_sig_otus, drug_label){
  
  #    drug <- "control"
  #    drug_sig_otus <- ""
  #    drug_label <- "No antibiotics"
  drug_rabund <- rel_abund_sig[metadata$abx == drug,]
  drug_metadata <- metadata[metadata$abx == drug,]
  
  n <- nrow(drug_rabund)
  
  drug_med <- apply(drug_rabund, 2, median)
  drug_uci <- apply(drug_rabund, 2, function(x){quantile(x, prob=0.75)})
  drug_lci <- apply(drug_rabund, 2, function(x){quantile(x, prob=0.25)})
  
  z <- barplot(drug_med, names.arg=rep("", length(drug_med)),
               ylim=c(0,1+max(drug_uci)), xlim=c(0,x_max), axes=F,
               col="white")
  
  warn_orig <- options("warn")$warn
  options("warn"= -1)
  
  arrows(x0=z, y0=drug_med, y1=drug_uci, angle=90, length=0.05)
  arrows(x0=z, y0=drug_med, y1=drug_lci, angle=90, length=0.05)
  
  options("warn"= warn_orig)
  
  text(x=z[sig_otus[o] %in% drug_sig_otus], y=-0.05*max(drug_uci),
       labels="*", cex=2, xpd=TRUE)
  
  axis(2, las=1)
  box()
  
  if(grepl(")", drug_label)){
    drug_label <- gsub(")", paste0("; N=", n, ")"), drug_label)
  } else {
    drug_label <- paste0(drug_label, "; N=", n, ")")
  }
  
  text(x=par("usr")[1], y=par("usr")[4]*1.175, label=drug_label,
       adj=c(0,1), cex=1.2, font=2, xpd=TRUE)
  
  
  
  
  #    summary_stats <- format(quantile(drug_metadata$CFU,
  #                            prob=c(0.25, 0.50, 0.75)), scientific=T, digits=2)
  #
  #    summary_string <- paste0(summary_stats[2], "~(", summary_stats[1], "-", summary_stats[3], ")")
  #
  #    summary_string <- gsub("(\\d\\.\\d*)e\\+0", "plain('\\1x10')^", summary_string)
  #    summary_string <- gsub("0e\\+00", "plain('<1x10')^2", summary_string)
  #
  #    summary_string <- paste0(summary_string, "~plain(' N=", n, "')")
  #
  #    text(x=par("usr")[2], y=1.05*par("usr")[4], labels=parse(text=summary_string),
  #                                adj=c(1,0), pos=2, cex=0.8, xpd=TRUE)
  
  summary_stats <- quantile(drug_metadata$CFU, prob=c(0.25, 0.50, 0.75))
  z2 <- barplot(summary_stats[2]+1, width=0.3, xlim=c(-0.1,0.5), ylim=c(1,1e9), log="y", axes=FALSE, col="white", names.arg="")
  
  warn_orig <- options("warn")$warn
  options("warn"= -1)
  
  arrows(x0=z2, y0=summary_stats[2]+1, y1=summary_stats[3]+1, angle=90, length=0.05)
  arrows(x0=z2, y0=summary_stats[2]+1, y1=summary_stats[1]+1, angle=90, length=0.05)
  
  options("warn"= warn_orig)
  
  box()
  axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))
  
  if(summary_stats[2] == 0){
    text(x=0.2, y=10, label=expression(plain('<10')^2))
  }
  
  z
  
}

#basically copied the rest of this from Pat and Alyx's repo

tiff(file="results/figures/figure1.tiff", width=4.5, height=10.0, unit="in", res=300)
par(cex=1.2)

layout_matrix <-matrix(c(
  1,  2,
  3,  4,
  5,  6,
  7,  8,
  9, 10,
  11, 12,
  13, 14,
  15, 16,
  17, 18
),nrow=9, byrow=T)
layout_matrix <- cbind(c(rep(19,8),0), layout_matrix, c(rep(20,8),0))

layout(layout_matrix, width=c(0.2, 1.25, 0.15, 0.2), height=c(rep(1,8),1.5))

#    par(mar=c(0.5,5,1.5,0.5))
par(mar=c(0.75,0.25,1.5,0.5), oma=c(0,0,0,0))

#deviates slightly here 
z <- single_drug_bars("none", "", "No antibiotics")

#not working
z <- single_drug_bars("amp", amp_sig_otus, "Ampicillin (0.5 mg/mL)")







