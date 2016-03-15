#Analysis starting 3 14 16

#lets first try Shannon diversity 
shannon <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.groups.summary", header = T)
plot(shannon$group, shannon$shannon)

#bring in metadata file, merge and subset to get some more interesting 
metadata <- read.table(file = "data/raw/abx_cdiff_metadata.tsv", header = T)
fullShannon <- merge(metadata, shannon)
#subset data by antibiotic
clindaShannon <- subset(fullShannon, abx == 'clinda')
#this looks like crap, do something different 
barplot(clindaShannon$shannon)

#ok instead lets try looking at simpson diversity. lets look at avg of cages per day over time per abx
simpson <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.groups.summary", header = T)
simpgroup <- simpson[, colnames(simpson) %in% c("group", "invsimpson")]
fullSimpson <- merge(metadata, simpgroup)

plot(blankplot)
for(i in 1:5){
day0 <- mean(fullSimpson[fullSimpson$cage==i&fullSimpson$day==0,'invsimpson']
day1 <- fullSimpson[fullSimpson$cage==i&fullSimpson$day==0,'invsimpson']

points(simp, type='l')
}
