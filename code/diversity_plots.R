## diversity vs colonization
#
##alpha diversity
#communities colonized to higher levels have lower diversity (alpha)
#	association between cfu and alpha (invsimpson and shannon - NS)
#	cfu vs # of otus (NS)
#	shared otus?
##beta diversity
#highly infected communities are most different than untreated
#	separation between untreated mice and all the highly infected communities (>1e6)
#communities that recover/elimnate cdifficile are more diverse
#	difference in diversity between highly infected 
#more change w/low diversity?
#more change with high cfu?
#
#need to remove dependence of daily sampling?



library(tidyverse)

metadata <- read.table(file = "data/raw/abx_cdiff_metadata.tsv", header = T)
alpha_df <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.ave-std.summary", header = T)

CFU_vs <- function(variable_name){
	data_frame <- metadata %>% 
		filter(day > 0, cdiff == TRUE) %>% 
		#filter(CFU > 0) %>% 
		select(group, CFU) %>% 
		mutate (CFU = CFU + 0.000000001) %>% 
		left_join(select(alpha_df, group, get(variable_name)))
	plot <- data_frame %>% 
		ggplot(aes(x = get(variable_name), y = CFU)) + 
			geom_point() + 
			scale_y_log10() + 
			geom_smooth(method = 'lm') +
			theme_bw() + 
			labs(x = variable_name, y = 'C.difficile CFU')
	fit_data <- summary(lm(get(variable_name) ~ CFU, data_frame))
	return(list(plot, fit_data))
}

CFU_vs('invsimpson')
CFU_vs('sobs')
CFU_vs('shannon')


#shared otus - among high or low?
1 get otus 
2 select samples that are in high (>1e5) and low (0)
3 

otus_df <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared", header = T)
tax_df <- read.table(file = "data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header = T)

#convert abundance to presence (T/F)
otu_presence <- otus_df %>%
	gather(otu, abundance, contains('Otu00')) %>%
	mutate(presence = ifelse(abundance > 0, T, F)) %>%
	select(Group, otu, presence) %>%
	spread(otu, presence)
# merge otu presence with colonization
high_low_df <- metadata %>%
	filter(day > 0, cdiff == TRUE) %>% 
	select(group, CFU) %>% rename(Group = group) %>% 
	filter(!between(CFU, 1, 1e5)) %>%
	mutate(colonization = ifelse(CFU < 1, F, T)) %>%
	inner_join(otu_presence) 
# get list of otus that are present in at least half of either group of samples
otu_list <- high_low_df %>%
	group_by(colonization) %>%
	summarise_each(funs(mean)) %>%
	ungroup() %>%
	select(contains('Otu00')) %>%
	summarise_each(funs(max)) %>%
	gather(otu, mean_presence) %>%
	filter(mean_presence >= 0.5) %>%
	select(otu)

cfu_otu_presence_df <- 
high_low_df %>%
	select(colonization, one_of(c(otu_list$otu))) %>%
	gather(otu, present, contains('Otu')) %>% 
	ggplot(aes(x = otu, y = count, group = colonization)) + 
		geom_point()



head(otus_df)
