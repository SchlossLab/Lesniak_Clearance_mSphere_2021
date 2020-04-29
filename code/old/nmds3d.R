library(ggplot2)
library(rgl)
library(magick)

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
nmds_3$abx <- as.character(nmds_3$abx)
nmds_3$dayplot <- nmds_3$day

nmds_3$dayplot[nmds_3$day < 0] <- -1

nmds_3$color <- as.numeric(factor(as.character(nmds_3$abx)))
nmds_3$size <- as.numeric(round(c(log10(nmds_3$CFU +1) + 1), 0))

# plot nmds
# plot in 3D
day_list <- split(nmds_3, nmds_3$dayplot)
abx_list <- split(nmds_3, nmds_3$abx)
cfu_df <- nmds_3[nmds_3$day > 0, ]
par3d(windowRect = c(0, 0, 700, 700))
with(nmds_3, 
	points3d(axis1, axis2, axis3, 
		#xlab = 'axis_1', ylab = 'axis_2', zlab = 'axis_3',
		width = 1000, size=NULL))
for(sz in c(1:10)){
	with(cfu_df[cfu_df$size == sz, ], 
			points3d(axis1, axis2, axis3, 
				#xlab = 'axis_1', ylab = 'axis_2', zlab = 'axis_3',
				width = 1000, size = sz)
				)
		#writeOBJ(paste0('test_abx_', unique(nmds_3$abx)[i], '.obj'), 
		#	#separateObjects = T
		#	pointRadius = c(seq(10, 40, 3.33)/1000)[sz]
		#	)
}

bgplot3d({
  plot.new()
  title(main = 'Antibiotic-treated conventional mice challenged with C. difficile', line = 3)
  mtext(side = 1, 'Points are sized by colonization level',
  	 line = 4, cex = 2)
  legend("topright", cex = 2,
  	pch=c(16,16,16,16), pt.cex = c(0.6,0.9,1.2,1.5), inset=c(0.02),
  	legend = c(expression(paste("10"^"2")), expression(paste("10"^"4")),
  		expression(paste("10"^"6")), expression(paste("10"^"8")))
  	)
})

#writeOBJ('test_legend.obj')
movie3d(spin3d(axis=c(0,0,1), rpm=3), dir = 'scratch/nmds3d/', duration=20, fps=10, movie="nmds3d_plot_by_cfu")
writeWebGL(filename = 'scratch/nmds3d/nmds_3d_by_cfu.html')