# File to generate validation data set
# Based on method and script used by 
# Kenta Suzuki, et al. 2017 Methods in Ecology and Evolution

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

# Generate Interaction matrix
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
a1 <- 7.5; a2 <- 2.7; a3 <- 2; B <- 0.1
GLVE <- function(interaction_matrix, K, gamma, noise_level, length, delta_time){
	time_series <- c()
	initial_state <- runif(numberofSpecies, min = 5, max = 9)
	current_state <- initial_state
	time_series <- c(0, current_state)
	for(tick in 1:length){
		new_state <- c()
		for(i in 1:numberofSpecies){
			xi <- current_state[i]
			if(gamma == 0){ 
				new_state[i] <- xi + (  
					( ri * (1 - sum(exp(current_state))/K) ) +  
					(sum(a1 * interaction_matrix[i, ] * exp(current_state))/ K ) ) * delta_time +  
					rnorm(1, mean = 0, sd = noise_level) 
			} else if(gamma == 1){ 
				new_state[i] <- xi + ( ( ri * (1 - sum(exp(current_state))/K) ) +  
					sum((a2 * interaction_matrix[i, ] * exp(current_state)^gamma)/ 
						(beta*K^gamma + exp(current_state)^gamma) ) ) * delta_time +  
					rnorm(1, mean = 0, sd = noise_level) 
			} else if(gamma >= 2){ 
				new_state[i] <- xi + ( ( ri * (1 - sum(exp(current_state))/K) ) +  
					sum((a3 * interaction_matrix[i, ] * exp(current_state)^gamma)/ 
						(beta*K^gamma + exp(current_state)^gamma) ) ) * delta_time +  
					rnorm(1, mean = 0, sd = noise_level) 
			} else {
				print('Error: gamma input incorrect')
			}
		}
		if(min(new_state) < 0) stop('Error: Species went extinct')
		time_series <- rbind(time_series, c(tick, new_state))
		current_state <- new_state
		}
	colnames(time_series) <- c('ticks', paste0('OTU_', 1:numberofSpecies))
	return(time_series)
}

interaction_matrix <- IM(numberofSpecies, connectance, cii, cij)	
time_series <- GLVE(interaction_matrix, K, gamma, noise_level, simulation_length, delta_time)
time_series %>%
	data.frame %>% 
	filter(ticks %in% seq((simulation_length - sampling_length), simulation_length, sampling_interval)) %>% 
	mutate(ticks = (ticks - (simulation_length - sampling_length))/sampling_interval) %>% 
	gather(OTU, abundance, contains('OTU')) %>%
	mutate(abundance = exp(abundance)) %>%  
	ggplot(aes(x = ticks, y = abundance, color = OTU)) + 
		geom_line()

# Create Sample Time Series



# Save data
write.table(time_series, paste0(save_dir, treatment_subset, '/validation_data_set.txt'), 
	quote = F, row.names = F)
