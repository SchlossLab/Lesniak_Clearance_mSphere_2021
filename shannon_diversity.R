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


for(treatment in levels(fullSimpson$abx)){
  plot(0, type='n', ylim=c(0,20), xlim=c(-6,10), xlab='day', ylab='inverse simpson', main=treatment)
  colors <- c('red','blue','green','pink','purple', 'orange')
  cages <- unique(fullSimpson[fullSimpson$abx==treatment, 'cage'])
  for(j in 1:length(cages)){
    days <- sort(unique(fullSimpson[fullSimpson$abx==treatment & fullSimpson$cage==cages[j],'day']))
    avg_simp <- c()
    for(i in 1:length(days)){
      avg_simp[i] <- mean(fullSimpson[fullSimpson$abx==treatment & fullSimpson$cage==cages[j] & fullSimpson$day==days[i], 'invsimpson'])
    }
    points(days, avg_simp, type='b', col=colors[j])
  }
  legend("top", legend=cages, col=colors, pch=1, lty=1, horiz=T)
} 
    
    
