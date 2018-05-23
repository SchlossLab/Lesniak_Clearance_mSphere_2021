library(multispatialCCM)
library(tidyverse)
library(cowplot)

#		test workflow with short simulated RPS data
#		setup to work with diff abx combinations
#		check into ientification of otus with 0 abundance

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F)
shared_file <- read.table(shared_file, sep = '\t', header = T)

#test with one abx treatment

#create vector of mouse by day stitched per abx (separate each mouse with NA)
#  begin with NA
Accm<-ccm_data_out$Accm
Bccm<-ccm_data_out$Bccm

#Calculate optimal E
maxE<-9 #Maximum E to test
#Matrix for storing output
Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")
#Loop over potential E values and calculate predictive ability
#of each process for its own dynamics
for(E in 2:maxE) {
#Uses defaults of looking forward one prediction step (predstep)
#And using time lag intervals of one time step (tau)
Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
}
#Look at plots to find E for each process at which
#predictive ability rho is maximized
matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
	xlab="E", ylab="rho", lwd=2)
legend("bottomleft", c("A", "B"), lty=1:2, col=1:2, lwd=2, bty="n")
#Results will vary depending on simulation.
# set maximum E for A and B
E_A<-2
E_B<-3

#Check data for nonlinear signal that is not dominated by noise
#Checks whether predictive ability of processes declines with
#increasing time distance
signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
	predsteplist=1:10)
signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
	predsteplist=1:10)

#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)
# Does A "cause" B?
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=100)
# Does B "cause" A?
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=100)

#Test for significant causal signal
CCM_significance_test<-ccmtest(CCM_boot_A,
	CCM_boot_B)

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot "A causes B"
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
	xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
	xlab="L", ylab="rho")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
	cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
		CCM_boot_A$rho+CCM_boot_A$sdevrho),
	lty=3, col=1)
#Plot "B causes A"
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
	cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
		CCM_boot_B$rho+CCM_boot_B$sdevrho),
	lty=3, col=2)
legend("topleft",
	c("A causes B", "B causes A"),
	lty=c(1,2), col=c(1,2), lwd=2, bty="n")


test_df <- meta_file %>% 
	filter(abx == 'clinda', cdiff == T, day >= 0) %>%
	select(group, cage, mouse, day, CFU) %>%
	left_join(select(shared_file, -label, -numOtus),
		by = c('group' = "Group")) %>%
	arrange(cage, mouse, day)%>%
	select(day, CFU, contains('Otu')) %>% 
	select(which(apply(test_df, 2, max) > 1 ))
NA_list <- which(test_df$day == 0)
test_df[test_df == 0] <-  sample(100,sum(test_df == 0, na.rm=T), replace = T)/100
test_df[NA_list, ] <- NA

set.seed(1)
output <- c()
for(i in 3:ncol(test_df)){
	Accm<-test_df$CFU
	Bccm<-test_df[,i]
	#Maximum E to test - one less than number of observations per sample
	maxE<-length(unique(test_df$day)) - 2 # one less for separating NAs and one less sample
	#Matrix for storing output
	Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")
	#Loop over potential E values and calculate predictive ability
	#of each process for its own dynamics
	for(E in 2:maxE) {
	#Uses defaults of looking forward one prediction step (predstep)
	#And using time lag intervals of one time step (tau)
	Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
	}
	#maximum E 
	E_A<-c(2:maxE)[which(Emat[,1] == max(Emat[,1], na.rm =T))]
	E_B<-c(2:maxE)[which(Emat[,2] == max(Emat[,2], na.rm =T))]
	embedding_dim_plot <- data.frame(cbind(Emat, E = c(2:maxE))) %>% 
		gather(OTU, rho, -E) %>% 
		left_join(data.frame(OTU = c('A', 'B'), Selected_E = c(E_A, E_B))) %>% 
		ggplot(aes(x = E, y = rho, color = OTU)) + 
			geom_line() + 
			geom_vline(aes(xintercept = Selected_E, color = OTU), linetype = 'dashed', size = 0.5) +
			labs(x = 'E', y = 'Pearson correlation coefficient (rho)', title = 'Best embedding dimension selection')
	#Check data for nonlinear signal that is not dominated by noise
	#Checks whether predictive ability of processes declines with
	#increasing time distance
	#See manuscript and R code for details
	signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
	predsteplist=1:10)
	signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
	predsteplist=1:10)
	prediction_step_plot <- rbind(data.frame(signal_A_out$predatout, OTU = 'A'),
		data.frame(signal_B_out$predatout, OTU = 'B')) %>% 
		ggplot(aes(x = predstep, y = rho, color = OTU)) + 
			geom_line() + 
			labs(x = 'Prediction Steps', y = 'Pearson correlation coefficient (rho)', 
				title = 'Predictive Power')
	#Run the CCM test
	#E_A and E_B are the embedding dimensions for A and B.
	#tau is the length of time steps used (default is 1)
	#iterations is the number of bootsrap iterations (default 100)
	# Does A "cause" B?
	#Note - increase iterations to 100 for consistant results
	CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=1000)
	# Does B "cause" A?
	CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=1000)
	#Test for significant causal signal
	#See R function for details
	CCM_significance_test<-ccmtest(CCM_boot_A, CCM_boot_B)
	current_ccm <- data.frame(t(CCM_significance_test), 
		cdiff_cause_otu = max(CCM_boot_A$rho), 
		otu_cause_cdiff = max(CCM_boot_B$rho),
		otu = colnames(test_df)[i],
		E_A = E_A,
		E_B = E_B)

	if(current_ccm$pval_b_cause_a <= 0.05){
		causal_otu <- as.character(current_ccm$otu)
		CCM_plot <- rbind(data.frame(causal = paste0('Cdiff_causes_', causal_otu), 
			lobs = CCM_boot_A$Lobs,
			rho = CCM_boot_A$rho,
			stdev_min = CCM_boot_A$rho - CCM_boot_A$sdevrho,
			stdev_max = CCM_boot_A$rho + CCM_boot_A$sdevrho),
		data.frame(causal = paste0(causal_otu, '_causes_Cdiff'), 
			lobs = CCM_boot_B$Lobs,
			rho = CCM_boot_B$rho,
			stdev_min = CCM_boot_B$rho - CCM_boot_B$sdevrho,
			stdev_max = CCM_boot_B$rho + CCM_boot_B$sdevrho )) %>% 
		#gather(level, value, rho, stdev_min, stdev_max) %>% 
		ggplot(aes(x = lobs)) + 
			geom_line(aes(y = rho, color = causal)) + 
			geom_ribbon(aes(ymin = stdev_min, ymax = stdev_max, fill = causal), alpha = 0.3) + 
			labs(x = 'L', y = 'Pearson correlation coefficient (rho)', color = '', fill = '') + 
			theme_bw() + 
			theme(legend.position="top", legend.direction="horizontal")

		dynamics_plot <- meta_file %>% 
			filter(abx == 'clinda', cdiff == T, day >= 0) %>%
			left_join(select(shared_file, -label, -numOtus),
				by = c('group' = "Group")) %>%
			select(cage, mouse, day, CFU, one_of(causal_otu)) %>% 
			gather(bacteria, counts, CFU, one_of(causal_otu)) %>% 
				ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), group = interaction(cage, mouse))) + 
					geom_line() + 
					facet_grid(bacteria~., scales = 'free_y') +
					theme_bw() + 
					labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', 
						title = 'Temporal Dynamics', subtitle = 'Colored by mouse') + 
					scale_x_continuous(breaks=seq(0,10, 1)) + 
					theme(legend.position = 'none')
			

		ggsave(paste0('scratch/ccm/ccm_cdiff_caused_by_', causal_otu, '.jpg'),
			plot_grid(embedding_dim_plot, prediction_step_plot, CCM_plot, dynamics_plot))

	}
	output <- rbind(output, current_ccm)
}

summary(output)

causal_otus <- output %>% 
	filter(pval_b_cause_a <= 0.05) %>% 
	arrange(desc(otu_cause_cdiff)) %>% 
	pull(otu) %>% 
	as.character

plot_causal_otus <- function(causal_otu){
	Accm <- test_df$CFU
	Bccm <- pull(test_df, causal_otu)
	E_A <- pull(output[output$otu == causal_otu, ], 'E_A')
	E_B <- pull(output[output$otu == causal_otu, ], 'E_B')
	#Run the CCM test
	# Does A "cause" B?
	CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=100)
	# Does B "cause" A?
	CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=100)

	CCM_plot <- rbind(data.frame(causal = paste0('Cdiff_causes_', causal_otu), 
			lobs = CCM_boot_A$Lobs,
			rho = CCM_boot_A$rho,
			stdev_min = CCM_boot_A$rho - CCM_boot_A$sdevrho,
			stdev_max = CCM_boot_A$rho + CCM_boot_A$sdevrho),
		data.frame(causal = paste0(causal_otu, '_causes_Cdiff'), 
			lobs = CCM_boot_B$Lobs,
			rho = CCM_boot_B$rho,
			stdev_min = CCM_boot_B$rho - CCM_boot_B$sdevrho,
			stdev_max = CCM_boot_B$rho + CCM_boot_B$sdevrho )) %>% 
		#gather(level, value, rho, stdev_min, stdev_max) %>% 
		ggplot(aes(x = lobs)) + 
			geom_line(aes(y = rho, color = causal)) + 
			geom_ribbon(aes(ymin = stdev_min, ymax = stdev_max, fill = causal), alpha = 0.3) + 
			labs(x = 'L', y = 'rho', color = '', fill = '') + 
			theme_bw() + 
			theme(legend.position="top", legend.direction="horizontal")

	dynamics_plot <- meta_file %>% 
		filter(abx == 'clinda', cdiff == T, day >= 0) %>%
		left_join(select(shared_file, -label, -numOtus),
			by = c('group' = "Group")) %>%
		select(cage, mouse, day, CFU, one_of(causal_otu)) %>% 
		gather(bacteria, counts, CFU, one_of(causal_otu)) %>% 
			ggplot(aes(x = day, y = counts, color = interaction(as.factor(mouse), as.factor(cage)), group = interaction(cage, mouse))) + 
				geom_line() + 
				facet_grid(bacteria~., scales = 'free_y') +
				theme_bw() + 
				labs(x = 'Day', y = 'Abundance \n (C difficle = CFU, Otu = 16s counts)', subtitle = 'Colored by mouse') + 
				scale_x_continuous(breaks=seq(0,10, 1)) + 
				theme(legend.position = 'none')
		

	ggsave(paste0('scratch/ccm/ccm_cdiff_caused_by_', causal_otu, '.jpg'),
		plot_grid(CCM_plot, dynamics_plot))
}

lapply(causal_otus, plot_causal_otus)