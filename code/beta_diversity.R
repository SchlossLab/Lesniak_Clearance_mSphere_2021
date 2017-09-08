library(ggplot2)
library(rgl)

cfu_labels <- read.table('data/mothur/CFU.design', header = T, sep = '\t')
cdiff_labels <- read.table('data/mothur/cdiff.design', header = T, sep = '\t')
abx_labels <- read.table('data/mothur/abx.design', header = T, sep = '\t')
cage_labels <- read.table('data/mothur/cage.design', header = T, sep = '\t')
day_labels <- read.table('data/mothur/day.design', header = T, sep = '\t')

# read in nmds axes file produced by mothur
time_nmds_3m <- read.table('data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.thetayc.0.03.lt.ave.nmds.3.axes', sep = '\t', header = T)

#merge nmds with metadata
nmds_3 <- merge(time_nmds_3m, 
	merge(cfu_labels,
		merge(abx_labels,
			merge(cdiff_labels,
				merge(cage_labels, day_labels)))))
nmds_3$CFU[is.na(nmds_3$CFU)] <- 0
nmds_3$abx[nmds_3$abx == 'vanc '] <- 'vanc'

# plot nmds
# plot 3 axis in 2D
ggplot(nmds_3, aes(x = axis1, y = axis2)) +
	geom_point(aes(color=as.factor(abx), size = CFU, shape = cdiff)) +
	theme_bw()
# plot in 3D
plot3d(nmds_3$axis1, nmds_3$axis2, nmds_3$axis3, col = factor(as.character(nmds_3$abx))

