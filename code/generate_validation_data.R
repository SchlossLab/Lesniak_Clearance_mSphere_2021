# File to generate validation data set
# Based on method and script used by 
# Kenta Suzuki, et al. 2017 Methods in Ecology and Evolution
#  https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.12814&file=mee312814-sup-0002-SupInfo.nb

library(tidyverse) # for tidyr/dplyr/%>% functions
library(deSolve) # ode() to run GLVE
library(purrr) # mp_dfr() to create replicate time series

input <- commandArgs(TRUE)
if(length(input) != 2) stop('Incomplete input. Please enter a seed and gamma (0, 1, or 2))')
seed <- as.numeric(input[1])
gamma <- as.numeric(input[2])
if(!is.numeric(seed)) stop('Please input a seed')
if(!is.numeric(gamma)) stop('Please input a gamma (0, 1, 2)')
set.seed(seed)

# Set simulation parameters
delta_time <- 1
simulation_length <- 1000
sampling_length <- 11
sampling_interval <- 1
replicates <- 12
subsample_level <- 2000

# Set model parameters
numberofSpecies <- 7
connectance <- 0.5 # proportion of realized interactions among potential ones
cii <- -0.5 #intraspecific interaction
cij_min <- 0.05 # interspecific interaction
cij_max <- 1 # interspecific interaction
K <- 10^4 # carrying capacity (upper limit of abundance, set to subsample amount)
#gamma <- 2 # select type of Holling Type Functional Response
# 0 for linear interactions, 1 for intermediate nonlinearity, >1 all cases nonlinearity
beta <- (0.1*K)^gamma # half-saturation constant of the interspecific interaction
sigma <- 0.1 # used 0.01, 0.1, 0.2
noise_level <-  sigma * (delta_time^0.5) # 
ri <- 1 # intrinsic growth rate
a1 <- 7.5; a2 <- 2.7; a3 <- 2#; B <- 0.1 #alternative to beta

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

# Subsample simulated count data
subsample <- function (count_data){
	subsample_df <- apply(select(count_data, -time), 
			1, function(x){
				values <- round(exp(x)) %>% 
					rep(names(.), .) %>% 
					sample(., subsample_level, replace = T)}) %>% 
		t %>% 
		data.frame(stringsAsFactors = F) %>% 
		mutate(time = count_data$time - min(count_data$time)) %>% 
		gather(draw, otu, -time) %>% 
		group_by(time) %>% 
		count(otu) %>% 
		spread(otu, n) %>% 
		ungroup
	subsample_df[is.na(subsample_df)] <- 0
	return(subsample_df)
}

date()
print(paste0('Using gamma = ', gamma, ' and seed set to ', seed))

# Create Sample Time Series
initial_state <- runif(numberofSpecies, min = 5, max = 9)
# Generate interaction matrix
# interaction_matrix <- as.matrix(
#	read.table('data/process/validation/validation_interaction_matrix_gamma2_seed1.txt'))
interaction_matrix <- IM(numberofSpecies, connectance, cii, cij)
# test and create matrix until interactions produce persistent existence
x <- NULL
while(!is.numeric(x)){
	x <- try(map_df(1:20, function(x){
			ode(y = runif(numberofSpecies, min = 5, max = 9), times = 0:5000, 
				func = GLVE, parms = list(noise_level = 0), method = "iteration") %>% 
				data.frame
			}), 
		silent = T)
	if (class(x)=="try-error") {
		#cat("ERROR1: ", x, "\n")
		Sys.sleep(1)
		#print("Trying new matrix")
		interaction_matrix <- IM(numberofSpecies, connectance, cii, cij)
	} else break 
}
print('Interaction Matrix set')
date()

# Generate simulation of time series
time_series <- map_dfr(1:10, function(x){
	# simulate time series
	ts <- ode(y = initial_state, # initial values (vector)
		times = 0:simulation_length, # time sequence desired
		func = GLVE, # R-function with func <- function(t, y, parms, ...), 
					# t = current time, y = current estimate, parms = parameters
					# must output a list of derivatives of y with respect to time
		parms = list(noise_level = noise_level), # list of parameters input to model
		method = "iteration") %>% 
		data.frame %>% 
		mutate(time = time * delta_time) %>% 
		# Use tail timepoints (after temporal dynamics have become stable)
		filter(time >= c((simulation_length - sampling_length * sampling_interval) * delta_time)) 
	# subsample time series
	ts <- subsample(ts) %>% 
		mutate(replicate = x)
	return(ts)
	})

colnames(time_series) <- gsub('X', 'OTU_', colnames(time_series))

suffix <- paste0('_gamma', gamma, '_seed',  seed, '.') 
ggsave(paste0('data/process/validation/validation_time_series', suffix, 'jpg'),
	time_series %>% 
		gather(OTU, abundance, contains('OTU')) %>%
		mutate(abundance = abundance) %>%  
		ggplot(aes(x = time, y = abundance, color = OTU, group = interaction(OTU, replicate))) + 
			geom_line())

# Save data
write.table(time_series, paste0('data/process/validation/validation_temporal_data', suffix, 'txt'),
	quote = F, row.names = F)
write.table(interaction_matrix, paste0('data/process/validation/validation_interaction_matrix', suffix, 'txt'),
	quote = F, row.names = F, col.names = F)
