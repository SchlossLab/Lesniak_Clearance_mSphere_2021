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


#shared otus - among high or no colonization?
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
#subset metadata to samples
presence_df <- metadata %>%
	filter(day > 0, cdiff == TRUE) %>% 
	select(group, CFU, day) %>% rename(Group = group) %>% 
	filter(!between(CFU, 1, 1e5)) %>%
	mutate(colonization = ifelse(CFU < 1, F, T)) %>%
	inner_join(otu_presence)
# get list of otus that are present in at least half of either group of samples
otu_list <- presence_df %>%
	select(colonization, contains('Otu00')) %>%
	gather(otu, presence, contains('Otu00'))%>%
	group_by(colonization, otu) %>%
	summarise_each(funs(mean)) %>%
	ungroup()  %>% 
	group_by(otu) %>%
	summarise_each(funs(max)) %>%
	filter(presence >= 0.5) %>%
	select(otu)
# merge otu presence with colonization
presence_by_day_df <- presence_df %>%
	select(colonization, day, one_of(otu_list$otu)) %>%
	gather(otu, presence, contains('Otu00')) %>%
	group_by(colonization, otu, day) %>%
	summarise_each(funs(mean)) %>%
	ungroup()
presence_by_day_df$otu_labels <- gsub("Otu0*", "", presence_by_day_df$otu)

#test for significant differences between colonized/uncolonized by otu
sig_presence <- c()
for(i in unique(presence_by_day_df$otu)){
	colonized <- presence_by_day_df %>%
		filter(otu == i & colonization == T) %>%
		select(presence)
	uncolonized <- presence_by_day_df %>%
		filter(otu == i & colonization == F) %>%
		select(presence)
	sig_presence <- rbind(sig_presence,
		list(i, as.numeric(wilcox.test(c(colonized$presence), c(uncolonized$presence))$p.value)))
}
significant_otu_list <- unlist(sig_presence[unlist(sig_presence[,2]) < 0.05/46, 1])

presence_by_day_df %>%
	ggplot(aes(x = otu, y = presence, color = colonization)) + 
		stat_summary(fun.data=median_hilow, position=position_dodge(0.2)) + #, conf.int=.5) +
		labs(title = 'Comparison of OTUs present in Colonized/Uncolonized Samples',
			subtitle = '(Median and IQR shown, Only OTUs present in 50% of samples shown)',
			y = 'Portion of Samples with OTU Present', x = 'OTU') + 
		scale_color_manual(labels = c("Uncolonized (C. difficle CFU = 0)", "Colonized (C. difficle CFU >= 1e5)"), values = c("blue", "red")) +
		scale_x_discrete(labels=unique(presence_by_day_df$otu_labels)) + 
		theme_bw() + 
		theme(legend.justification=c(0,0), legend.position=c(0.005,0.005), 
			legend.title = element_blank(),   legend.box.background = element_rect(),
			legend.box.margin = margin(0, 0, 0, 0), legend.key.height = unit(0.3, "cm"),
			axis.text.x = element_text(angle = 45, hjust = 1))


# merge otu presence with colonization
all_presence_by_day_df <- presence_df %>%
	select(colonization, day, contains('Otu00')) %>%
	gather(otu, presence, contains('Otu00')) %>%
	group_by(colonization, otu, day) %>%
	summarise_each(funs(mean)) %>%
	ungroup()
all_presence_by_day_df$otu_labels <- gsub("Otu0*", "", all_presence_by_day_df$otu)

#test for significant differences between colonized/uncolonized by otu
all_sig_presence <- c()
for(i in unique(all_presence_by_day_df$otu)){
	colonized <- all_presence_by_day_df %>%
		filter(otu == i & colonization == T) %>%
		select(presence)
	uncolonized <- all_presence_by_day_df %>%
		filter(otu == i & colonization == F) %>%
		select(presence)
	all_sig_presence <- rbind(all_sig_presence,
		list(i, as.numeric(wilcox.test(c(colonized$presence), c(uncolonized$presence))$p.value)))
}
all_significant_otu_list <- unlist(all_sig_presence[unlist(all_sig_presence[,2]) < 0.05/nrow(all_sig_presence), 1])

all_sig_by_day_df <- all_presence_by_day_df %>%
	filter(otu %in% all_significant_otu_list) 

all_sig_by_day_df %>%
	ggplot(aes(x = otu, y = presence, color = colonization)) + 
		stat_summary(fun.data=median_hilow, position=position_dodge(0.2), fun.args=list(conf.int=0.5)) +
		labs(title = 'Comparison of OTUs present in Colonized/Uncolonized Samples',
			subtitle = '(Median and IQR shown, Only OTUs present in 50% of samples shown)\nAll OTUs plotted are significant (p < 0.05) after BF multiple comparisons correction ',
			y = 'Portion of Samples with OTU Present', x = 'OTU') + 
		scale_color_manual(labels = c("Uncolonized (C. difficle CFU = 0)", "Colonized (C. difficle CFU >= 1e5)"), values = c("blue", "red")) +
		scale_x_discrete(labels=unique(all_sig_by_day_df$otu_labels)) + 
		theme_bw() + 
		theme(legend.justification=c(0,0), legend.position=c(0.005,0.005), 
			legend.title = element_blank(),   legend.box.background = element_rect(),
			legend.box.margin = margin(0, 0, 0, 0), legend.key.height = unit(0.3, "cm"),
			axis.text.x = element_text(angle = 45, hjust = 1))


