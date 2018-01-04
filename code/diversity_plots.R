# diversity vs colonization

#alpha diversity
communities colonized to higher levels have lower diversity (alpha)
	association between cfu and alpha
	cfu vs # of otus
	shared otus?
#beta diversity
highly infected communities are most different than untreated
	separation between untreated mice and all the highly infected communities (>1e6)
communities that recover/elimnate cdifficile are more diverse
	difference in diversity between highly infected 
more change w/low diversity?
more change with high cfu?

need to remove dependence of daily sampling?



library(tidyverse)

metadata <- read.table(file = "data/raw/abx_cdiff_metadata.tsv", header = T)
alpha_df <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.ave-std.summary", header = T)

invsim_df <- metadata %>% 
	filter(CFU > 0) %>% 
	select(group, CFU) %>% 
	left_join(select(alpha_df, group, invsimpson))

invsim_df %>% 
	ggplot(aes(x = invsimpson, y = CFU)) + 
		geom_point() + 
		scale_y_log10() + 
		geom_smooth(method = 'lm')

summary(lm(invsimpson ~ CFU, invsim_df))