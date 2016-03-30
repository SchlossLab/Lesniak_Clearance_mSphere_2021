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



#Or the ddply way. turn this into a script/loop this weekend, dont forget to add SD lines
#test <- ddply(fullSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))

library(plyr)
plot(0, type='n', ylim=c(0,25), xlim=c(0,10), xlab='day', ylab='inverse simpson', main='Diversity over time, all treatments')

clindaSimpson <- subset(fullSimpson, abx == 'clinda')
strepSimpson <- subset(fullSimpson, abx == 'strep')
ampSimpson <- subset(fullSimpson, abx == 'amp')
cefSimpson <- subset(fullSimpson, abx == 'cef')
cipSimpson <- subset(fullSimpson, abx == 'cipro')
metSimpson <- subset(fullSimpson, abx == 'metro')
vancSimpson <- subset(fullSimpson, abx == 'vanc')

avg_clin_simp <- ddply(clindaSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))
avg_strep_simp <- ddply(strepSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))
avg_amp_simp <- ddply(ampSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))
avg_cef_simp <- ddply(cefSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))
avg_cip_simp <- ddply(cipSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))
avg_met_simp <- ddply(metSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))
avg_vanc_simp <- ddply(vancSimpson,~day,summarise, avg = mean(invsimpson), sd = sd(invsimpson))

points(avg_strep_simp$day, avg_strep_simp$avg, type = 'b',lwd = 2,  col = 'red')
points(avg_clin_simp$day, avg_clin_simp$avg, type = 'b',lwd = 2,  col = 'orange')
points(avg_amp_simp$day, avg_amp_simp$avg, type = 'b', lwd = 2, col = 'lightblue')
points(avg_cef_simp$day, avg_cef_simp$avg, type = 'b', lwd = 2, col = 'green')
points(avg_cip_simp$day, avg_cip_simp$avg, type = 'b', lwd = 2, col = 'blue')
points(avg_met_simp$day, avg_met_simp$avg, type = 'b', lwd = 2, col = 'purple')
points(avg_vanc_simp$day, avg_vanc_simp$avg, type = 'b', lwd = 2, col = 'brown')

legend("top", legend=c('strep', 'clinda', 'amp', 'cef', 'cipro', 'metro', 'vanc'), col=c('red', 'orange', 'yellow', 'green', 'blue', 'purple', 'brown'), cex=0.7, pch=1, lty=1, lwd = 2, horiz=T)

arrows(avg_vanc_simp$day, avg_vanc_simp$avg-avg_vanc_simp$sd, avg_vanc_simp$day, avg_vanc_simp$avg+avg_vanc_simp$sd, length=0.05, angle=90, code=3, col = 'brown')
arrows(avg_strep_simp$day, avg_strep_simp$avg-avg_strep_simp$sd, avg_strep_simp$day, avg_strep_simp$avg+avg_strep_simp$sd, length=0.05, angle=90, code=3, col = 'red')
arrows(avg_clin_simp$day, avg_clin_simp$avg-avg_clin_simp$sd, avg_clin_simp$day, avg_clin_simp$avg+avg_clin_simp$sd, length=0.05, angle=90, code=3, col = 'orange')
arrows(avg_amp_simp$day, avg_amp_simp$avg-avg_amp_simp$sd, avg_amp_simp$day, avg_amp_simp$avg+avg_amp_simp$sd, length=0.05, angle=90, code=3, col = 'lightblue')
arrows(avg_cef_simp$day, avg_cef_simp$avg-avg_cef_simp$sd, avg_cef_simp$day, avg_cef_simp$avg+avg_cef_simp$sd, length=0.05, angle=90, code=3, col = 'green')
arrows(avg_cip_simp$day, avg_cip_simp$avg-avg_cip_simp$sd, avg_cip_simp$day, avg_cip_simp$avg+avg_cip_simp$sd, length=0.05, angle=90, code=3, col = 'blue')
arrows(avg_met_simp$day, avg_met_simp$avg-avg_met_simp$sd, avg_met_simp$day, avg_met_simp$avg+avg_met_simp$sd, length=0.05, angle=90, code=3, col = 'purple')


#now to make the OTU table do this..  but maybe Pat doesn't want this rightnow 
otu_table <- read.table(file = "data/mothur/genus.1.subsample.shared.trimmed", header = T)