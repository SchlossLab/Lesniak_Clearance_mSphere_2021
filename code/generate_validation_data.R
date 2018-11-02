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
noise <- 0.0316
interaction <- matrix(unlist(read.table('~/True_interaction_matrix.txt')),
	byrow = T, nrow = 10)

# Pat ODE
library(deSolve)

ode_model <- function(t, x, p){
	n <- length(x)
	data <- ( p$growth*(1-sum(x)/1) + p$interaction %*% x ) * 
			(1 + rnorm(n, mean = 0, sd = p$noise))
	list(data * x)
}

generate_community <- function(initial, time, growth, interaction=NULL,
														perturbation=NULL, susceptibility=NULL,
														n_replicates=1, noise=0, method="ode45"){
	output_data <- list()
	if(is.list(initial)){
		n_replicates <- length(initial)
	}
	for(replicate in 1:n_replicates){
		params <- list(growth = growth, interaction = interaction,
											noise = noise)
		start <- NULL
		if(is.list(initial)){
			start <- initial[[replicate]]
		} else {
			start <- initial
		}
		data <- ode(start, time, ode_model, params, method)[,-1]
		output_data[[replicate]] <- as.matrix(data)
	}
	names(output_data) <- 1:n_replicates
	output_data
}

growth <- rep(0.5,numberofSpecies)
tmp <- runif(numberofSpecies)
initial <- tmp/sum(tmp)
time <- seq(1,1000,10)
times <- list(time, time, time, time, time, time, time, time, time, time)
names(times) <- as.character(length(times))
X_obs <- generate_community(initial=initial, time=time, growth=growth,
							interaction=interaction, n_replicates=length(times), noise=noise)
X_obs$`1` %>% data.frame %>% mutate(ticks = time) %>% gather(OTU, abundance, -ticks) %>% mutate(OTU = gsub('X', 'OTU_', OTU)) %>% ggplot(aes(x = ticks, y = abundance, color = OTU)) + geom_line()
# issue: limited dynamics, not affected by noise, replicates don't differ

# Suzuki 2017 (https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.12814&file=mee312814-sup-0002-SupInfo.nb)

GLVE <- function(interaction_matrix, K, gamma, noise_level, length, delta_time){
	time_series <- c()
	initial_state <- runif(numberofSpecies, min = 5, max = 9)
	current_state <- initial_state
	time_series <- c(0, current_state)
	beta <- (0.1*K)^gamma # half-saturation constant of the interspecific interaction
	#NOISE <- rnorm(1, mean = 0, sd = noise_level) 
	for(tick in 1:10){
		new_state <- c()
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
#		for(i in 1:numberofSpecies){
#			xi <- current_state[i]
#			if(gamma == 0){ 
#				new_state[i] <- xi + (  
#					( ri * (1 - sum(exp(current_state))/K) ) +  
#					(sum(a1 * interaction_matrix[i, ] * exp(current_state))/ K ) ) * delta_time +  
#					rnorm(1, mean = 0, sd = noise_level) 
#			} else if(gamma == 1){ 
#				new_state[i] <- xi + ( ( ri * (1 - sum(exp(current_state))/K) ) +  
#					sum((a2 * interaction_matrix[i, ] * exp(current_state)^gamma)/ 
#						(beta + exp(current_state)^gamma) ) ) * delta_time +  
#					rnorm(1, mean = 0, sd = noise_level) 
#			} else if(gamma >= 2){ 
#				new_state[i] <- xi + ( ( ri * (1 - sum(exp(current_state))/K) ) +  
#					sum((a3 * interaction_matrix[i, ] * exp(current_state)^gamma)/ 
#						(beta + exp(current_state)^gamma) ) ) * delta_time +  
#					rnorm(1, mean = 0, sd = noise_level) 
			} else {
				print('Error: gamma input incorrect')
			}
#		}
#		current_state - new_state
		if(min(new_state) < 0) stop('Error: Species went extinct')
		time_series <- rbind(time_series, c(tick, new_state))
		current_state <- new_state
		}
	colnames(time_series) <- c('ticks', paste0('OTU_', 1:numberofSpecies))
	return(time_series)
}
numberofSpecies <- 7
interaction_matrix <- IM(numberofSpecies, connectance, cii, cij)	
interaction_matrix <- matrix(unlist(read.table('~/ssm_interaction_matrix.txt')),
	byrow = T, nrow = 7)
time_series <- GLVE(interaction_matrix, K, gamma, 0, simulation_length, delta_time)
time_series %>%
	data.frame %>% 
	filter(ticks %in% seq((simulation_length - sampling_length), simulation_length, sampling_interval)) %>% 
	mutate(ticks = (ticks - (simulation_length - sampling_length))/sampling_interval) %>% 
	gather(OTU, abundance, contains('OTU')) %>%
	mutate(abundance = exp(abundance)) %>%  
	ggplot(aes(x = ticks, y = abundance, color = OTU)) + 
		geom_line()
# without noise, eventually becomes stable, with noise, continually interact


# Fisher Mehta - LIMITS (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4108331/bin/pone.0102451.s001.nb)
interaction <- matrix(unlist(read.table('~/True_interaction_matrix.txt')),
	byrow = T, nrow = 10)
numberofSpecies <- 10
L <- 100
random_time_series <- function(interaction_matrix, sigma, L){
	# simulates a timeseries of length L for interaction matrix interaction_matrix and 
	# stochasticity sigma
	tc <- matrix(nrow = c(L + 1), ncol = c(numberofSpecies))
	x <- matrix(rlnorm((L + 1) * numberofSpecies, mean = 0, sd = sigma), nrow = c(L + 1))
	tc[1,] <- x[1,]
	for(i in 1:L){
		tc[i + 1, ] <- x[c(i + 1), ] * tc[i, ] * exp(interaction %*% tc[i, ])
	}
	return(tc)
}
LIMITS_ts <- random_time_series(interaction, noise_level, L)
apply(LIMITS_ts, 1, function(x)x/sum(x)) %>% t %>% 
	data.frame %>% 
	mutate(ticks = 0:L - 35) %>% 
	gather(OTU, abundance, contains('X')) %>% 
	mutate(OTU = gsub('X','OTU_', OTU),		
		abundance = abundance) %>% 
	filter(ticks >= 0) %>% 
	ggplot(aes(x = ticks, y = abundance, color = OTU)) + 
		geom_line()

# Create Sample Time Series



# Save data
write.table(time_series, paste0(save_dir, treatment_subset, '/validation_data_set.txt'), 
	quote = F, row.names = F)