#---------------------------------------------------------------------#
# 4. Creating model ####
#
# This creates the model, simple as that.
# 
# 
#---------------------------------------------------------------------#

# Import raw data
raw_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/species_model_data.RDS"))

# Just a couple of transformation functions. Greta, and most of these Bayesian packages, requires standardised data. I've log-transformed two of the variables as well, SCI and buffer_5000m_population_2006.
log_sc_gr <- function(vector) {new_vector <- as.numeric(scale(log(vector+1)))
return(new_vector)}
sc_gr <- function(vector) {new_vector <- as.numeric(scale(vector))
return(new_vector)}

# Define which distance threshold we're using
raw_data$no_n_pop <- raw_data[,paste0("no_n_pop_",population_threshold,"km")]


# We now get rid of all lakes in the native range and transform variables.
env_data <- raw_data %>%
  filter(native==0) %>%
  transmute(lake_area = log_sc_gr(area_km2),
            distance_to_road = log_sc_gr(distance_to_road),
            temperature = sc_gr(eurolst_bio10),
            distance_nearest_population = log_sc_gr(dist_n_pop),
            human_footprint = sc_gr(HFP),
            number_nearby_populations = log_sc_gr(no_n_pop),
            upstream_presence = sc_gr(upstream_presence),
            downstream_presence = sc_gr(downstream_presence))

# Get rid of any parameters we want to ignore.
if(!is.na(parameters_to_ignore)) {
  env_data <- env_data[,-parameters_to_ignore]
}

# Create presence data
intro_data <- raw_data %>%
  filter(native==0) %>%
  dplyr::select(introduced)

id_data <- raw_data %>%
  filter(native==0) %>%
  dplyr::select(locationID)


# Turn them into Greta arrays
env_dataGr <- as_data(env_data)
intro_dataGr <- as_data(intro_data)

# A couple of values that will make things easier.
n_env <- ncol(env_dataGr)
n_sites <- nrow(env_dataGr)
n_int <- length(interaction_terms)



print("Data is downloaded, parameters are ready")

# Create priors
alpha <- normal(0, 10)

beta <- normal(0, 10, dim=n_env+n_int)

# Prepare interaction terms
if (!is.na(interaction_terms)) {
  for (i in 1:length(interaction_terms)) {
  first_set <- interaction_terms[[i]]
  first_term <- env_dataGr[,first_set[1]]
  second_term <- env_dataGr[,first_set[2]]
  first_x_second <- first_term * second_term
  if (i == 1) {
    interaction_term <- first_x_second %*% beta[n_env+i,]
  } else {interaction_term <- interaction_term + first_x_second %*% beta[n_env+i,]}
}
}

# Defines our equation
if (!is.na(interaction_terms)) {
  linear_predictor <- alpha + env_dataGr %*% beta[1:n_env,] +
    interaction_term
} else {
  linear_predictor <- alpha + env_dataGr %*% beta
  }


# These two give our equation the required transformations
p <- ilogit(linear_predictor)
distribution(intro_dataGr) <- bernoulli(p)

# These define and visualise our model.
prelim_model <- model(beta)
# plot(prelim_model)

print("Constructed model, running draws now.")

# The following then creates our MCMC draws
whole_draws <- greta::mcmc(prelim_model,n_samples = 500, warmup = 500)
print("Finished running first draws, running extra draws now.")
whole_draws_extra <- extra_samples(whole_draws,n_samples = 500)

print("Finished running draws, saving data now.")

whole_model_output <- list(draws = whole_draws_extra, beta = beta, alpha = alpha, p = p)
saveRDS(whole_model_output, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                        "/whole_model_output_",population_threshold,"km.RDS"))

whole_model_data <- list(raw_data = raw_data, env_data = env_data, intro_data = intro_data, id_data = id_data)
saveRDS(whole_model_data, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                      "/whole_model_data_",population_threshold,"km.RDS"))
