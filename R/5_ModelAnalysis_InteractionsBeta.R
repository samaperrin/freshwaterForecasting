
#---------------------------------------------------------------------#
# 5. Run Model Analysis
#---------------------------------------------------------------------#


# This will actually be a script for checking either model

model_output <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_output_",population_threshold,"km.RDS"))
model_data_extra <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_data_",population_threshold,"km.RDS"))

# Define which distance threshold we're using
raw_data <- model_data_extra$raw_data
raw_data$no_n_pop <- raw_data[,paste0("no_n_pop_",population_threshold,"km")]

draws <- model_output$draws
beta <- model_output$beta
p <- model_output$p
alpha <- model_output$alpha



#---------------------------------------------------------------------#
# 5A. Run Model Convergence Diagnostics
#---------------------------------------------------------------------#

### At this point let's run some quick convergence diagnostics. Values for gelman 
#   diagnostics should be close to 1, for the effective size the higher the better.
gelmans <- coda::gelman.diag(calculate(beta,values = draws),multivariate = FALSE)
number_converged <- nrow(gelmans$psrf[gelmans$psrf[,1] < 1.1,])
print(paste0(number_converged," of ",nrow(gelmans$psrf)," parameters converged."))

effective_sizes <- summary(coda::effectiveSize(calculate(beta,values=draws)))


#---------------------------------------------------------------------#
# 5B. Get Model Deviance
#---------------------------------------------------------------------#

# Whatever data/model we are using for out parameters, we are now applying the parameters to all the data.
# That means we need to scale all the data using whatever means and SDs with which we scaled the data that
# was used for the moel.
data_4scaling <- raw_data %>% filter(native == 0)
#env_var <- c("HFP", "nearby_pops", "pop_dist", "area_km2", "perimeter_m", "distance_to_road", "eurolst_bio10")
data_4scaling <- data_4scaling %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temp = eurolst_bio10,
            pop_distL = log(dist_n_pop+1),
            HFP = HFP,
            n_pop = log(no_n_pop+1),
            upstream_presence = upstream_presence,
            downstream_presence = downstream_presence)

if(!is.na(parameters_to_ignore)) {
  data_4scaling <- data_4scaling[,-parameters_to_ignore]
}

# Get the means and sds to work with
means_relData <- apply(data_4scaling, 2, mean)
SDs_relData <- apply(data_4scaling, 2, sd)

env_data_takeMeans <- t(apply(data_4scaling, 1, "-", means_relData))
env_data_scaled <- t(apply(env_data_takeMeans, 1, "/", SDs_relData))

if (!is.na(interaction_terms[[1]][1])) {
  for (i in 1:length(interaction_terms)) {
    first_set <- interaction_terms[[i]]
    first_term <- as.data.frame(env_data_scaled[,first_set[1]])
    second_term <- as.data.frame(env_data_scaled[,first_set[2]])
    if (i == 1) {interaction_term <- first_term * second_term
    colnames(interaction_term) <- paste0(first_set[1],'x',first_set[2])
    } else {
      additional_term <- first_term * second_term
      colnames(additional_term) <- paste0(first_set[1],'x',first_set[2])
      interaction_term <- cbind(interaction_term, additional_term)
    }
  }
  env_data_scaled <- cbind(env_data_scaled, interaction_term)
}

print("Data is scaled.")

### Get betas and alphas
calc_beta <- apply(as.matrix(calculate(beta,values=draws)),2,mean)
calc_alpha <- apply(as.matrix(calculate(alpha,values=draws)),2,mean)

# Get interaction terms
# Prepare interaction terms
if (!is.na(interaction_terms[[1]][1])) {
  for (i in 1:length(interaction_terms)) {
    first_set <- interaction_terms[[i]]
    first_term <- env_data_scaled[,first_set[1]]
    second_term <- env_data_scaled[,first_set[2]]
    first_x_second <- first_term * second_term
    if (i == 1) {
      interaction_term <- first_x_second * calc_beta[length(means_relData)+i]
    } else {interaction_term <- interaction_term + first_x_second * calc_beta[length(means_relData)+i]}
  }
}


# Get probabilities of colonisation for lakes we have introductions for
attempt <- as.matrix(env_data_scaled) %*% as.matrix(calc_beta)
  
eta <- sweep(attempt, 2, calc_alpha, "+")
expeta <- exp(eta)
init_probabilities <- expeta/(1+expeta)


p_dev <- init_probabilities[,1]
Y_dev <- model_data_extra$intro_data$introduced
model_deviance <- deviance_yp(Y_dev,p_dev)

print("Calculated deviance. Now going through getting uncertainty.")


#---------------------------------------------------------------------#
# 5C. Get Uncertainty Values
#---------------------------------------------------------------------#

# Create the matrix to stuff all the values into
intVal_mat <- data.frame(lower=numeric(),
                         mean=numeric(), 
                         upper=numeric(), 
                         stringsAsFactors=FALSE) 

if (nrow(env_data_scaled) < 10000) {
  int_eta <- pred_env(env_data_scaled,alpha,beta)
  int_val <- ilogit(int_eta)
  int_draws <- calculate(int_val,values=draws)
  
  # The following creates our probabilities of introduction.
  int_pred_ints <- as.data.frame(t(apply(as.matrix(int_draws) , 2 , quantile , probs = c(0.025,0.5,0.975) , na.rm = TRUE )))
  intVal_mat <- rbind(intVal_mat,int_pred_ints)
} else {
  loops <- ceiling(nrow(env_data_scaled)/10000)
  
  time1 <- Sys.time()
  for (i in 1:(loops)) {
    if(i != (loops)) {
      analyse_df_chunk <- env_data_scaled[(1+(i-1)*10000):(i*10000),]
    } else {
      analyse_df_chunk <- env_data_scaled[(1+(i-1)*10000):nrow(env_data_scaled),]
    }
    int_eta <- pred_env(analyse_df_chunk,alpha,beta)
    int_val <- ilogit(int_eta)
    int_draws <- calculate(int_val,values=draws)
    
    # The following creates our probabilities of introduction.
    int_pred_ints <- as.data.frame(t(apply(as.matrix(int_draws) , 2 , quantile , probs = c(0.025,0.5,0.975) , na.rm = TRUE )))
    intVal_mat <- rbind(intVal_mat,int_pred_ints)
    
    # Quick function to let us know how the loop is progressing
    time2 <- Sys.time()
    difftime_1 <- round(as.numeric(difftime(time2, time1,
                                            units = "mins")),4)
    if (i %% 10 == 0) {print(paste0("Run ", i, " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*(loops+1)/i-difftime_1,10), " minutes left.") )}
  }
  
}


colnames(intVal_mat) <- c("lower","mean","upper")

intVal_mat$width <- with(intVal_mat,(upper-lower)/2)

#---------------------------------------------------------------------#
# 5D. Get Parameter Estimates
#---------------------------------------------------------------------#

## Now let's check the beta estimates
beta_ints <- get_beta_list(draws,beta_shared=beta,species_names="introduced",
                           env_names=colnames(env_data_scaled))


model_analysis <- list(intervals = intVal_mat, deviance = model_deviance, 
                       conv_diags = list(gelmans, effective_sizes), parameter_estimates = beta_ints)

saveRDS(model_analysis,paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                              "/whole_model_analytics_",population_threshold,"km_interaction.RDS"))

