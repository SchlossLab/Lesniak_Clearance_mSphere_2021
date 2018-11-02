# File to generate validation data set
# Based on method and script used by 
# Kenta Suzuki, et al. 2017 Methods in Ecology and Evolution
#  https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.12814&file=mee312814-sup-0002-SupInfo.nb

library(tidyverse)
library(deSolve)

seed <- input
set.seed(seed)

# Set simulation parameters
delta_time <- 0.1
simulation_length <- 5000
sampling_length <- 1000
sampling_interval <- 10

# Set model parameters
numberofSpecies <- 7
connectance <- 0.5 # proportion of realized interactions among potential ones
cii <- -0.5 #intraspecific interaction
cij_min <- 0.05 # interspecific interaction
cij_max <- 1 # interspecific interaction
K <- 10^4 # carrying capacity (upper limit of abundance, set to subsample amount)
gamma <- 2 # select type of Holling Type Functional Response
# 0 for linear interactions, 1 for intermediate nonlinearity, >1 all cases nonlinearity
beta <- (0.1*K)^gamma # half-saturation constant of the interspecific interaction
sigma <- 0.1 # used 0.01, 0.1, 0.2
noise_level <-  sigma * (delta_time^0.5) # 
ri <- 1 # intrinsic growth rate
a1 <- 7.5; a2 <- 2.7; a3 <- 2#; B <- 0.1

# Create function to generate Interaction matrix
IM <- function(numberofSpecies, connectance, cii, cij){ 
	# create a vector containing random sample of interactions
	total_interactions <- numberofSpecies^2
	interactions <- runif(total_interactions, cij_min, cij_max) * sample(c(-1,1), total_interactions, replace = T)
	# reduce interaction to level of connectance
	interactions[sample(length(interactions), (numberofSpecies^2 - numberofSpecies)*(connectance))] <- 0
	mat <- matrix(interactions,	nrow = numberofSpecies, ncol = numberofSpecies)
	diag(mat) <- cii
	return(mat)
}

# Create Generalized Lotka-Voltera Equation Function
initial_state <- runif(numberofSpecies, min = 5, max = 9)

GLVE <- function(t, current_state, p){ 
	noise_level <- p$noise_level
	if(gamma == 0){
	new_state <- current_state + ( 
			( ri * (1 - sum(exp(current_state))/K) ) + 
			( interaction_matrix %*% ((a1 * exp(current_state))/ K)) ) * delta_time + 
			rnorm(numberofSpecies, mean = 0, sd = noise_level)
	} else if(gamma == 1){
		new_state <- current_state + ( 
			( ri * (1 - sum(exp(current_state))/K) ) + 
			( interaction_matrix %*% c((a2 * exp(current_state)^gamma)/ 
				(beta + exp(current_state)^gamma)) )) * delta_time + 
			rnorm(numberofSpecies, mean = 0, sd = noise_level)
	} else if(gamma >= 2){
		new_state <- current_state + ( 
			( ri * (1 - sum(exp(current_state))/K) ) + 
			( interaction_matrix %*% c((a3 * exp(current_state)^gamma)/ 
				(beta + exp(current_state)^gamma)) )) * delta_time + 
			rnorm(numberofSpecies, mean = 0, sd = noise_level)
	} else {
		print('Error: gamma input incorrect')
	}
	if(min(new_state) < 0) stop('Error: Species went extinct')
	list(new_state)
}

# Create Sample Time Series
interaction_matrix <- IM(numberofSpecies, connectance, cii, cij)	
#interaction <- matrix(unlist(read.table('~/True_interaction_matrix.txt')),
#	byrow = T, nrow = 10)

params <- list()

time_series <- ode(y = initial_state, # initial values (vector)
	times = 0:10, # time sequence desired
	func = GLVE, # R-function with func <- function(t, y, parms, ...), 
				# t = current time, y = current estimate, parms = parameters
				# must output a list of derivatives of y with respect to time
	parms = list(noise_level = 0), # list of parameters input to model
     method = "iteration")

time_series <- GLVE(interaction_matrix, K, gamma, 0, 10, delta_time)
time_series %>%
	data.frame %>% 
	#filter(ticks %in% seq((simulation_length - sampling_length), simulation_length, sampling_interval)) %>% 
	#mutate(ticks = (ticks - (simulation_length - sampling_length))/sampling_interval) %>% 
	gather(OTU, abundance, -time) %>%
	mutate(ticks = time) %>% 
	mutate(abundance = exp(abundance)) %>%  
	ggplot(aes(x = ticks, y = abundance, color = OTU)) + 
		geom_line()

# Save data
write.table(time_series, paste0(save_dir, treatment_subset, '/validation_temporal_data.txt'), 
	quote = F, row.names = F)
write.table(interaction_matrix, paste0(save_dir, treatment_subset, '/validation_interaction_matrix.txt'), 
	quote = F, row.names = F)