#---------------------------------------------------------------------#
# 4. Creating model ####
#
# This creates the model, simple as that.
# 
# 
#---------------------------------------------------------------------#

# We are going to run validation data, right from the start

k <- 5

standardising_values <- readRDS("./Data/standardising_values.RDS")

means_locs <- standardising_values$means
sds_locs <- standardising_values$SDs

# First thing to do is set up all the data

data_for_species <- list()

for (i in 1:length(species_list)) {
  raw_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),"/species_model_data.RDS"))
  
  if (use_weighted_absences == TRUE) {
    data_for_splitting <- raw_data %>%
      filter(native==0) %>%
      filter(weighted_absences == 1 | introduced == 1)
  } else {
  data_for_splitting <- raw_data %>%
    filter(native==0)
  }
  
  # Now that we have all the data that will be used, we can split the data
  presences <- split(which(data_for_splitting$introduced == 1),1:k)
  absences <- split(which(data_for_splitting$introduced == 0),1:k)
  
  #Now let's standardise all the data
  data_for_splitting$no_n_pop <- data_for_splitting[,paste0("no_n_pop_",population_threshold,"km")]
  standardised_data <- data_for_splitting %>%
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
  scaled_data_almost <- sweep(standardised_data[,2:10],2,means_locs,"-")
  scaled_data_all <-cbind(standardised_data[,c(1,12:14)],  sweep(scaled_data_almost,2,sds_locs,"/"))
  
  # Now split everything up
  species_split_list <- list()
  for (r in 1:k) {
    validation_data <- scaled_data_all[c(presences[[r]],absences[[r]]),]
    training_data <- scaled_data_all[!(scaled_data_all$locationID %in% validation_data$locationID),]
    species_split_list[[r]] <- list(validation_data = validation_data, training_data = training_data)
  }
  data_for_species[[species_list[i]]] <- species_split_list
}

n_env_full <- 9 - length(parameters_to_ignore) + length(interaction_terms)

kfold_params <- list()

for (j in 1:k) {
  
  # first step is to split up data
  
  
  # Create priors
  global_alpha <- normal(0, 10, dim = 1)
  global_alpha_sd <- uniform(0, 10, dim = 1) 
  alpha <- normal(global_alpha, global_alpha_sd, dim = length(species_list))
  alpha_list <- cbind(alpha,alpha,alpha,alpha,alpha)
  
  global_betas <- normal(0, 10, dim = n_env_full)
  global_betas_sd <- uniform(0, 10, dim = n_env_full)
  beta <- cbind(normal(global_betas, global_betas_sd, dim = c(n_env_full,1)),
                normal(global_betas, global_betas_sd, dim = c(n_env_full,1)),
                normal(global_betas, global_betas_sd, dim = c(n_env_full,1)),
                normal(global_betas, global_betas_sd, dim = c(n_env_full,1)),
                normal(global_betas, global_betas_sd, dim = c(n_env_full,1)))
  beta_list <- list(beta,beta,beta,beta,beta)
  
  p_list <- list()
  
  # We are comparing the effects of variables across species. This is great, but seeing as each species has different
  # absences and presences as a result of being found in different lakes, when we scale the parameters they'll end up
  # having differenk <- 5t values unless we used a common measurement to scale them. Hence what we do now.
  for (s in 1:length(species_list)) {
    
    relevant_training_data <- data_for_species[[species_list[s]]][[1]]$training_data
    env_data_all <- relevant_training_data[,5:13]
    
    
    # Get rid of any parameters we want to ignore.
    if(!is.na(parameters_to_ignore)) {
      env_data <- env_data_all[,-parameters_to_ignore]
    } else {env_data <- env_data_all}
    
    # Create presence data
    intro_data <- relevant_training_data %>%
      dplyr::select(introduced)
    
    id_data <- relevant_training_data %>%
      dplyr::select(locationID)
    
    
    # Turn them into Greta arrays
    env_dataGr <- as_data(env_data)
    intro_dataGr <- as_data(intro_data)
    
    # A couple of values that will make things easier.
    n_env <- ncol(env_dataGr)
    n_env_interaction <- ncol(env_dataGr) + length(interaction_terms)
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
    p_list[[s]] <- ilogit(linear_predictor)
    distribution(intro_dataGr) <- bernoulli(p_list[[s]])
    
    
    
  }
  
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
  
  
  kfold_model_params <- list(draws = whole_draws_extra, beta = beta, alpha = alpha,
                             globals = list(global_alpha = global_alpha, global_betas = global_betas),
                             p = p_list)
  
  kfold_params[[j]] <- kfold_model_params
}

# These define and visualise our model.
# plot(prelim_model)

print("Constructed model, running draws now.")

  
print("Finished running draws, saving data now.")

# Figure out data file name

  
whole_model_output <- list(model_ouput = kfold_params, split_data = data_for_species)
if (use_weighted_absences == TRUE) {
  saveRDS(whole_model_output, file=paste0("./Data/whole_model_output_pseudoAbs_kfold_pooled.RDS"))
} else {
  saveRDS(whole_model_output, file=paste0("./Data/whole_model_output_kfold_pooled.RDS"))
}
  