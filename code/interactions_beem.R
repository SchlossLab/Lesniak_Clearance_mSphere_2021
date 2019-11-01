# Try to implement BEEM
# https://github.com/csb5/BEEM
# C Li, K R Chng, J S Kwah, T V Av-Shalom, L Tucker-Kellogg & N Nagarajan. (2019). An expectation-maximization algorithm enables accurate ecological modeling using longitudinal metagenome sequencing data. Microbiome.
#https://link.springer.com/epdf/10.1186/s40168-019-0729-z?shared_access_token=ie-5Gr1CKndobAUnYsRcf2_BpE1tBhCbnbw3BuzI2RPFn4dHsUppvFZnytFP_i4dN4fjxZldoW9nTbD9IBYoDiyG0aTOZ1-blvHrz2pi8MZjSVrljY30ecdNpBPboVg_0VEqW-f-AXv5HxR7njzpHxNCdP8z7T2gvmzjIiKdi3U%3D

library(beem)

# beem metadata
# sampleID    isIncluded    subjectID    measurementID
#    sampleID: sample IDs matching the first row of the OTU table
#    isIncluded: whether the sample should be included in the analysis (1-include, 0-exclude)
#    subjectID: indicator for which biological replicate the sample belongs to
#    measurementID: time in standardized units from the start of the experiment
meta <- meta_file %>%
	filter(day > 0, cdiff == T, abx == 'strep') %>% 
	mutate(subjectID = paste(cage, mouse, sep = '_'),
		measurementID = day + abs(min(day)),
		isIncluded = 1) %>% 
	left_join(sample_id, by = c('group')) %>% 
	select(group, sampleID, isIncluded, subjectID, measurementID) 

sample_id <- data.frame(group = meta$group, sampleID = 1:length(meta$group))

beem_counts <- shared_file %>% 
	rownames_to_column %>% 
	left_join(sample_id, by = c(rowname = 'group')) %>% 
	select(sampleID, contains('Otu')) %>% 
	t

# tutorial metadata
# meta <- read.table('~/BEEM/vignettes/props_et_al_analysis/metadata.sel.txt', head =T)
meta$purterbID <- 0
# tutorial counts
# df.counts <- read.table('~/BEEM/vignettes/props_et_al_analysis/counts.sel.txt', head =T,row.names=1, sep='\t')
df.counts <- beem_counts[-1, ]
colnames(df.counts) <- sample_id$sampleID
## tss counts
df.counts.tss <- apply(df.counts, 2, function(x) x/sum(x))
## keep high abundant otus
fil <- apply(df.counts.tss,1,mean)>0.001
df.counts.sel <- df.counts.tss[fil,]
## scale counts to a defined median
scaling <- 1e5
dat <- df.counts.sel/median(colSums(df.counts.sel)) * scaling

## Run BEEM
res <- EM(dat=dat, meta=meta)
### kept getting to the following error, but unclear where it is coming from 
# Preprocessing data ...
# Error in { : 
#   task 1 failed - "task 1 failed - "task 1 failed - "invalid 'times' argument"""

## Estimate parameters
biomass <- biomassFromEM(res)
write.table(biomass, 'biomass.txt', col.names=F, row.names=F, quote=F)
gLVparameters <- paramFromEM(res, counts, metadata)
write.table(gLVparameters, 'gLVparameters.txt', col.names=T, row.names=F, sep='\t' , quote=F)