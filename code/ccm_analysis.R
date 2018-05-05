library(multispatialCCM)
library(tidyverse)

meta_file   <- 'data/process/abx_cdiff_metadata_clean.txt'
shared_file <- 'data/mothur/abx_time.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared'
meta_file   <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F)
shared_file <- read.table(shared_file, sep = '\t', header = T)

test with one abx treatment
create vector of mouse by day stitched per abx 

test_df <- meta_file %>% 
	filter(abx == 'clinda', cdiff == T, day >= 0) %>%
	select(group, cage, mouse, day, CFU) %>%
	left_join(select(shared_file, -label, -numOtus),
		by = c('group' = "Group")) %>%
	arrange(cage, mouse, day)
test_df[test_df$day == 0, ] <- NA

for 
Accm<-test_df$CFU
Bccm<-test_df$Otu000012

#Simulate data to use for multispatial CCM test
#See function for details - A is causally forced by B,
#but the reverse is not true.
ccm_data_out<-make_ccm_data()
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
#Using the seed we provide,
#maximum E for A should be 2, and maximum E for B should be 3.
#For the analyses in the paper, we use E=2 for all simulations.
E_A<-2
E_B<-3
#Check data for nonlinear signal that is not dominated by noise
#Checks whether predictive ability of processes declines with
#increasing time distance
#See manuscript and R code for details
signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
predsteplist=1:10)
signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
predsteplist=1:10)
#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)
# Does A "cause" B?
#Note - increase iterations to 100 for consistant results
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=100)
# Does B "cause" A?
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=100)
#Test for significant causal signal
#See R function for details
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
	select(day, CFU, contains('Otu'))
NA_list <- which(test_df$day == 0)
test_df[test_df == 0] <-  sample(100,sum(test_df == 0, na.rm=T), replace = T)/100
test_df[NA_list, ] <- NA

output <- c()
for(i in 2:ncol(test_df)){
	Accm<-test_df$CFU
	Bccm<-test_df[,7]
	maxE<-6 #Maximum E to test
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
	#Check data for nonlinear signal that is not dominated by noise
	#Checks whether predictive ability of processes declines with
	#increasing time distance
	#See manuscript and R code for details
	signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
	predsteplist=1:10)
	signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
	predsteplist=1:10)
	#Run the CCM test
	#E_A and E_B are the embedding dimensions for A and B.
	#tau is the length of time steps used (default is 1)
	#iterations is the number of bootsrap iterations (default 100)
	# Does A "cause" B?
	#Note - increase iterations to 100 for consistant results
	CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=100)
	# Does B "cause" A?
	CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=100)
	#Test for significant causal signal
	#See R function for details
	CCM_significance_test<-ccmtest(CCM_boot_A, CCM_boot_B)
	output <- rbind(output, data.frame(t(CCM_significance_test), 
		cdiff_cause_otu = max(CCM_boot_A$rho), 
		otu_cause_cdiff = max(CCM_boot_B$rho),
		otu = colnames(test_df)[i]))
}

