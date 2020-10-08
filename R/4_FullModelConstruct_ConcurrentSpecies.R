#---------------------------------------------------------------------#
# 4. Creating model ####
#
# This creates the model, simple as that.
# 
# 
#---------------------------------------------------------------------#

# Import raw data

# Create priors
global_alpha <- normal(0, 10, dim = 1)
global_alpha_sd <- uniform(0, 10, dim = 1) 
alpha <- normal(global_alpha, global_alpha_sd, dim = length(species_list))

n_env <- 9 - length(parameters_to_ignore) + length(interaction_terms)

global_betas <- normal(0, 10, dim = n_env)
global_betas_sd <- uniform(0, 10, dim = n_env)
beta <- cbind(normal(global_betas, global_betas_sd, dim = c(n_env,1)),
              normal(global_betas, global_betas_sd, dim = c(n_env,1)),
              normal(global_betas, global_betas_sd, dim = c(n_env,1)),
              normal(global_betas, global_betas_sd, dim = c(n_env,1)),
              normal(global_betas, global_betas_sd, dim = c(n_env,1)))

p_list <- list()
species_data <- list()

# We are comparing the effects of variables across species. This is great, but seeing as each species has different
# absences and presences as a result of being found in different lakes, when we scale the parameters they'll end up
# having different values unless we used a common measurement to scale them. Hence what we do now.

species_data_import <- list()
raw_data_import <- list()
for (i in 1:length(species_list)) {
  raw_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),"/species_model_data.RDS"))
  raw_data_import[[i]] <- raw_data
  data_for_weighting <- raw_data %>%
    filter(native==0)
  data_for_weighting$no_n_pop <- data_for_weighting[,paste0("no_n_pop_",population_threshold,"km")]
  species_data_import[[i]] <- data_for_weighting %>%
    transmute(locationID = locationID,
              areaL = log(area_km2+1),
              dist_roadL = log(distance_to_road+1),
              temp = eurolst_bio10,
              pop_distL = log(dist_n_pop+1),
              HFP = HFP,
              n_pop = log(no_n_pop+1),
              upstream_presence = upstream_presence,
              downstream_presence = downstream_presence,
              nearby_weighted = log(nearby_weightings*10000+1),
              native, introduced, presence, weighted_absences)
}

# Need to use different means and SDs for different variables. Some variables are consistent across 
# species, but others (downstream presence, nearby extant populations) will vary from species to
# species.
standardising_values <- readRDS("./Data/standardising_values.RDS")

means_locs <- standardising_values$means
sds_locs <- standardising_values$SDs

#Now, start putting together species models one by one

for (s in 1:length(species_list)) {
  raw_data_inModel <- species_data_import[[s]]
  
  # Define which distance threshold we're using
  if (use_weighted_absences == TRUE) {
    data_for_model <- raw_data_inModel %>%
      dplyr::filter(introduced == 1 | weighted_absences == 1)
  } else {data_for_model <- raw_data_inModel}
  
  
  # We now get rid of all lakes in the native range and transform variables.
  env_data_all_almost <- sweep(data_for_model[,2:10],2,means_locs,"-")
  env_data_all <- sweep(env_data_all_almost,2,sds_locs,"/")
  
  # Get rid of any parameters we want to ignore.
  if(!is.na(parameters_to_ignore)) {
    env_data <- env_data_all[,-parameters_to_ignore]
  } else {env_data <- env_data_all}
  
  # Create presence data
  intro_data <- data_for_model %>%
    dplyr::select(introduced)
  
  id_data <- data_for_model %>%
    dplyr::select(locationID)
  
  
  # Turn them into Greta arrays
  env_dataGr <- as_data(env_data)
  intro_dataGr <- as_data(intro_data)
  
  # A couple of values that will make things easier.
  n_env <- ncol(env_dataGr)
  n_sites <- nrow(env_dataGr)
  n_int <- length(interaction_terms)
  
  print("Data is downloaded, parameters are ready")
  
  
  
  # Prepare interaction terms
  if (!is.na(interaction_terms)) {
    for (i in 1:length(interaction_terms)) {
      first_set <- interaction_terms[[i]]
      first_term <- env_data_all[,first_set[1]]
      second_term <- env_data_all[,first_set[2]]
      first_x_second <- first_term * second_term
      if (i == 1) {
        interaction_term <- first_x_second %*% beta[n_env+i,s]
      } else {interaction_term <- interaction_term + first_x_second %*% beta[n_env+i,s]}
    }
  }
  
  # Defines our equation
  if (!is.na(interaction_terms)) {
    linear_predictor <- alpha[s,] + env_dataGr %*% beta[1:n_env,s] +
      interaction_term
  } else {
    linear_predictor <- alpha[s,] + env_dataGr %*% beta[,s]
  }
  
  
  # These two give our equation the required transformations
  p_list[[s]] <- iprobit(linear_predictor)
  distribution(intro_dataGr) <- bernoulli(p_list[[s]])
  
  species_model_data <- list(raw_data = raw_data_inModel, env_data = env_data, intro_data = intro_data, id_data = id_data)
  
  species_data[[species_list[s]]] <- species_model_data
  
}

# These define and visualise our model.
prelim_model <- model(beta)
# plot(prelim_model)

print("Constructed model, running draws now.")

# The following then creates our MCMC draws

Lmin <- 10
Lmax <- round(Lmin*1.5)

whole_draws <- greta::mcmc(prelim_model,n_samples = initial_runs, warmup = initial_runs,
                           sampler = hmc(Lmin = Lmin, Lmax = Lmax))
print("Finished running first draws, running extra draws now.")
whole_draws_extra <- extra_samples(whole_draws,n_samples = n_extra_runs)

whole_model_params <- list(draws = whole_draws_extra, beta = beta, alpha = alpha,
                           globals = list(global_alpha = global_alpha, global_betas = global_betas),
                           p = p_list)


print("Finished running draws, saving data now.")

# Figure out data file name

if (use_weighted_absences == FALSE) {
  if (use_weighted_distances == FALSE) {
    file_distance_measure <- paste0(population_threshold,"km_fullAbs") 
  } else { 
    file_distance_measure <- "weightedDist_fullAbs" 
  }
} else {
  if (use_weighted_distances == FALSE) {
    file_distance_measure <- paste0(population_threshold,"km_pseudoAbs") 
  } else { 
    file_distance_measure <- "weightedDist_pseudoAbs" 
  }
}
  
whole_model_output <- list(species_data = species_data, whole_model_params = whole_model_params)
saveRDS(whole_model_output, file=paste0("./Data/whole_model_output_",file_distance_measure,"_pooled.RDS"))

