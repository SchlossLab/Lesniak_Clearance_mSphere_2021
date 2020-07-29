#Analysis starting 3 14 16

#lets first try Shannon diversity 
shannon <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.ave-std.summary", header = T)
plot(shannon$group, shannon$shannon)

#bring in metadata file, merge and subset to get some more interesting 
metadata <- read.table(file = "data/raw/abx_cdiff_metadata.tsv", header = T)
fullShannon <- merge(metadata, shannon)
#subset data by antibiotic
clindaShannon <- subset(fullShannon, abx == 'clinda')
#this looks like crap, do something different 
barplot(clindaShannon$shannon)

#ok instead lets try looking at simpson diversity. lets look at avg of cages per day over time per abx
simpson <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.ave-std.summary", header = T)
simpgroup <- simpson[, colnames(simpson) %in% c("group", "invsimpson")]
fullSimpson <- merge(metadata, simpgroup)

#this will make a diversity plot for each treatment that averages the mice in each cage and plots diversity per cage over time 
for(treatment in levels(fullSimpson$abx)){
  plot(0, type='n', ylim=c(0,20), xlim=c(-6,10), xlab='day', ylab='inverse simpson', main=treatment)
  colors <- c('red','blue','green','pink','purple', 'orange', 'brown', 'black', 'yellow')
  cages <- unique(fullSimpson[fullSimpson$abx==treatment, 'cage'])
  for(j in 1:length(cages)){
    days <- sort(unique(fullSimpson[fullSimpson$abx==treatment & fullSimpson$cage==cages[j],'day']))
    avg_simp <- c()
    for(i in 1:length(days)){
      avg_simp[i] <- mean(fullSimpson[fullSimpson$abx==treatment & fullSimpson$cage==cages[j] & fullSimpson$day==days[i], 'invsimpson'])
    }
      points(days, avg_simp, type='b', col=colors[j])
      
  }
  legend("top", legend=cages, col=colors, cex=0.7, pch=1, lty=1, horiz=T)
} 
 

#to make a plot that will avg all of the mice together in a treatment and plot over time lets do this   



#the ddply way. it's a loop now!!! 
library(plyr)
plot(0, type='n', ylim=c(0,25), xlim=c(0,10), xlab='day', ylab='inverse simpson', main='Diversity over time, all treatments')
new_colors <- c('red','orange','turquoise3','green','blue', 'purple', 'brown', 'black')
j <- 1
treatments <- c()
for(i in levels(fullSimpson$abx)){
  if(i != 'none'){
  temp <- subset(fullSimpson, abx == i)
  avg_i_simp <- ddply(temp, ~day, summarise, avg = mean(invsimpson), sd = sd(invsimpson))
  points(avg_i_simp$day, avg_i_simp$avg, type = 'b', lwd = 2, col = new_colors[j])
  arrows(avg_i_simp$day, avg_i_simp$avg-avg_i_simp$sd, avg_i_simp$day, avg_i_simp$avg+avg_i_simp$sd, length=0.05, angle=90, code=3, col = new_colors[j])
  treatments[j] <- i
  j <- j+1
  }
}
legend("top", legend=treatments, col=new_colors, cex=0.7, pch=1, lty=1, lwd = 2, horiz=T)


#now to make the OTU table do this..  but maybe Pat doesn't want this rightnow 
otu_table <- read.table(file = "data/mothur/genus.1.subsample.shared.trimmed", header = T)