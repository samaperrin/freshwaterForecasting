#---------------------------------------------------------------------#
# 6. Forecasting
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 6A. Get Scaling Parameters
#
# Variables within the model will change, but they have to be scaled 
# using the same mands and SDs as our original data, so we need to define
# those first.
#
#---------------------------------------------------------------------#


# Import data and draws and define components

model_output <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                               "/whole_model_output_",population_threshold,"km.RDS"))
model_data_extra <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                   "/whole_model_data_",population_threshold,"km.RDS"))
catchment_data <- readRDS("./Data/lakesInCatchments.RDS")
catchment_mapDF <- map_df(catchment_data$catchmentsByLake, ~data.frame(waterBodyID=.),.id="ebint")


raw_data <- model_data_extra$raw_data

draws <- model_output$draws
beta <- model_output$beta
p <- model_output$p
alpha <- model_output$alpha

# Whatever data/model we are using for out parameters, we are now applying the parameters to all the data.
# That means we need to scale all the data using whatever means and SDs with which we scaled the data that
# was used for the moel.
model_data <- model_data_extra$raw_data %>%
  filter(native == 0)
#env_var <- c("HFP", "nearby_pops", "pop_dist", "area_km2", "perimeter_m", "distance_to_road", "eurolst_bio10")
data_4scaling <- model_data %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temp = eurolst_bio10,
            pop_distL = log(dist_n_pop+1),
            HFP = HFP,
            n_pop = log(no_n_pop+1),
            upstream_presence = upstream_presence,
            downstream_presence = downstream_presence)

# Get the means and sds to work with
means_relData <- apply(data_4scaling, 2, mean)
SDs_relData <- apply(data_4scaling, 2, sd)

# Now we get our already scaled data
env_data <- model_data_extra$env_data

if (!is.na(interaction_terms)) {
  for (i in 1:length(interaction_terms)) {
    first_set <- interaction_terms[[i]]
    first_term <- env_data[,first_set[1]]
    second_term <- env_data[,first_set[2]]
    if (i == 1) {interaction_term <- first_term * second_term
    } else {
      interaction_term <- cbind(interaction_term, first_term * second_term)
    }
  }
  env_data <- cbind(env_data, interaction_term)
}

print("Data is scaled.")

#---------------------------------------------------------------------#
# 6B. Set Figures Beforehand
#
# We need to set up a bunch of parameters to use in the forecasting
# script. 
#
#---------------------------------------------------------------------#

### Get betas and alphas
calc_beta <- apply(as.matrix(calculate(beta,values=draws)),2,mean)
calc_alpha <- apply(as.matrix(calculate(alpha,values=draws)),2,mean)

# Get betas and alphas for the looping
beta_init <- as.matrix(calculate(beta,values=draws))
alpha_init <- as.matrix(calculate(alpha,values=draws))

set.seed(147)
chosen_iterations <- sample(1:nrow(alpha_init), n_loops)

beta_mat <- beta_init[chosen_iterations,]
alpha_mat <- alpha_init[chosen_iterations,]

# Need vector of where we have presences outside of the native range at the moment
presences <- model_data_extra$intro_data$introduced

nn_all <- get.knnx(raw_data[,c("utm_x","utm_y")],model_data[c("utm_x","utm_y")],k=100)
nn_all$nn.index[nn_all$nn.dist > population_threshold * 1000 | nn_all$nn.dist == 0] <- 0


periods <- list()

# These variables define actions in the loop
number_reps <- 5

# This defines in incremental increase over 50 years based on a 1 degree increase in temperature over 50 years.
# 10 represents a 1 degree increase. If the temp_scenariop is unnecessary, it just adds 0 every time.
temp_step <- ifelse(temp_scenario == TRUE,21/SDs_relData["temp"],rep(0,length(SDs_relData["temp"])))

# These variables will be used to scale projections of numnber of populations nearby and distance to nearest population
mean_n_pop <- means_relData["n_pop"]
sd_n_pop <- SDs_relData["n_pop"]

sd_dist <- SDs_relData["pop_distL"]
mean_dist <- means_relData["pop_distL"]

mean_upstream <- means_relData["upstream_presence"]
sd_upstream <- SDs_relData["upstream_presence"]

mean_downstream <- means_relData["downstream_presence"]
sd_downstream <- SDs_relData["downstream_presence"]

introductions <- matrix(NA,nrow=nrow(env_data),ncol=1)
periods <- matrix(NA,nrow=nrow(env_data),ncol=n_loops)
# introduction_probs <- list()
# populations <- list()

print("Other stuff in place, starting forecast run.")


#---------------------------------------------------------------------#
# 6C. Run Forecasting Loop
#
#
#---------------------------------------------------------------------#


time1 <- Sys.time()
# data <- list()
for (s in 1:n_loops) {
  for (j in 1:number_reps) {
    # So, first step is to create predictions based purely off a rise in temperature. Because this is slightly different to what we'll do in the other steps, there needs to be an if clause.
    
    if (j==1) {
      # We create a table containing our new data
      
      tempIncrease <- rnorm(nrow(env_data), temp_step/5, 0.05)
      newData <- env_data
      newData[,"temperature"] <- newData[,"temperature"] + tempIncrease
      
      # Update interaction terms
      if (!is.na(interaction_terms)) {
        for (i in 1:length(interaction_terms)) {
          first_set <- interaction_terms[[i]]
          first_term <- newData[,first_set[1]]
          second_term <- newData[,first_set[2]]
          newData[,length(means_relData)+i] <- first_term * second_term
        }
      }
      
      # Unfortunately because we're predicting for 600k + lakes, it's impossible to run the preditions all at one. So the loop below runs them in chunks of 10k lakes at a time.
      
      attempt <- as.matrix(newData) %*% as.matrix(beta_mat[s,])
      eta <- greta::sweep(attempt, 2, alpha_mat[s], "+")
      expeta <- exp(eta)
      probabilities <- expeta/(1+expeta)
      
      # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
      new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
      all_presences <- ifelse(presences ==1, 1, ifelse(new_presences == 1, 1, 0))
      
      newData2 <- newData
    } else {
      
      ### And now we move on to the second part of the loop, whereby we redefine those variables pertaining to closest population and numbers of close populations at each step
      
      # Introduce the new temperature thing first, it's the simplest
      tempIncrease <- rnorm(nrow(newData2), temp_step/5, 0.05)
      newData2[,"temperature"] <- newData2[,"temperature"] + tempIncrease
      
      # Now figure out the new distance to closest population
      pop_proximity <- cbind(all_presences, model_data[,c("utm_x","utm_y","decimalLatitude","decimalLongitude","locationID")])
      
      # Need to bind our new presences together with all the presences in the native range.
      data_presences_nonnative <- as.data.frame(pop_proximity %>% filter(all_presences == 1) %>% 
                                                  dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude))
      data_presences_native <- raw_data %>% filter(native == 1 & presence == 1) %>% 
                                                dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude)
      data_presences <- rbind(data_presences_nonnative, data_presences_native)
      
      data_all <- pop_proximity %>%
        dplyr::select(utm_x,utm_y,locationID, decimalLongitude, decimalLatitude) %>%
        distinct() %>%
        as.data.frame()
      
      # The get.knnx function returns the distance from a lake in table B to the k closest lakes in table A
      
      nn <- get.knnx(data_presences[c("utm_x","utm_y")],data_all[c("utm_x","utm_y")],k=2)
      
      # We then figure out the distance to the closest lake by taking the first column
      dist_to_closest_pop <- ifelse(nn$nn.dist[,1]==0,nn$nn.dist[,2],nn$nn.dist[,1])
      locationID <- data_all$locationID
      distance_data <- as.data.frame(cbind(dist_to_closest_pop,locationID))
      
      ordered_distances <- merge(model_data,distance_data,all.x=TRUE,by="locationID") %>%
        arrange(no_vatn_lnr)
      
      new_distances <- log(as.numeric(as.character(ordered_distances$dist_to_closest_pop))+1)
      
      # So now that we have the new measurements for closest population, these need to be scaled against those that we had for the initial population. So we subtract the mean and divide by the standard deviation.
      newData2[,"distance_nearest_population"] <- (new_distances-mean_dist)/sd_dist
      
      # Now we need to get the new measurements for number of close populations
      all_presence_table <- raw_data[,c("locationID","waterBodyID","no_vatn_lnr")]
      all_presence_table$present <- ifelse(all_presence_table$locationID %in% data_presences$locationID, 1, 0)
      
      rowNumbers_withPresence <- which(all_presence_table$present==1)
      nn_inds <- nn_all$nn.index
      nn_inds[!(nn_inds %in% rowNumbers_withPresence)] <- 0
      nn_inds[nn_inds %in% rowNumbers_withPresence] <- 1
      number_nearby_pop_vec <- log(rowSums(nn_inds)+1)
      # 
      # Scale it like we did for the last variable
      
      newData2[,"number_nearby_populations"] <- (number_nearby_pop_vec - mean_n_pop)/sd_n_pop
      
      # Now we need the number of populations upstream
      catchment_mapDF_ts <- catchment_mapDF
      
      # Get list of all water bodies with presence
      all_WB_present <- all_presence_table[all_presence_table$present == 1,"waterBodyID"]
      catchment_mapDF_ts$present <- ifelse(catchment_mapDF_ts$waterBodyID %in% all_WB_present, 1,0)
      # Define how many populations there are in each catchment
      new_upstream_pops <- catchment_mapDF_ts %>%
        group_by(ebint) %>%
        summarize(upstream_pops = sum(present))
      upstream_forNewData <- merge(model_data[,c("waterBodyID","ebint","no_vatn_lnr")], new_upstream_pops, all.x=TRUE, by = "ebint") %>%
        arrange(no_vatn_lnr)
      # We still have NAs where a lake isn't in its corresponding catchment
      # These can become zeroes, as there clearly isn't anything upstream in that catchment anyway
      upstream_forNewData$upstream_pops[is.na(upstream_forNewData$upstream_pops)] <- 0
      # If we have a situation where we have 2 presences and one of them is the focal lake, we need to reduce
      # the number of upstream lakes to 1.
      upstream_forNewData$upstream_pop_consol <- ifelse(upstream_forNewData$waterBodyID %in% all_WB_present,
                                                        upstream_forNewData$upstream_pops - 1, upstream_forNewData$upstream_pops)
      # We may have a couple of stray negatives. These can simply become 0.
      upstream_forNewData$upstream_pop_consol[upstream_forNewData$upstream_pop_consol<0] <- 0
      
      upstream_presence <- ifelse(upstream_forNewData$upstream_pop_consol == 0,0,1)
      newData2[,"upstream_presence"] <- (upstream_presence - mean_upstream)/sd_upstream
      
      # Need the number of populations downstream now
      basin_mapDF <- raw_data[,c("waterBodyID","eb_waterregionID")] 
      basin_mapDF <- basin_mapDF[!duplicated(basin_mapDF),]
      basin_mapDF$present <- ifelse(basin_mapDF$waterBodyID %in% all_WB_present, 1,0)
      # Define how many populations there are in each basin
      new_downstream_pops <- basin_mapDF %>%
        group_by(eb_waterregionID) %>%
        summarize(downstream_pops = sum(present))
      downstream_forNewData <- merge(model_data[,c("waterBodyID","eb_waterregionID","no_vatn_lnr")], new_downstream_pops
                                     , all.x=TRUE, by = "eb_waterregionID") %>%
        arrange(no_vatn_lnr)
      downstream_forNewData$downstream_pop_consol <- ifelse(downstream_forNewData$waterBodyID %in% all_WB_present,
                                                          downstream_forNewData$downstream_pops - 1, downstream_forNewData$downstream_pops)
      downstream_forNewData$downstream_pop_consol <- downstream_forNewData$downstream_pop_consol - upstream_forNewData$upstream_pop_consol
      # Once more, get rid of negatives
      downstream_forNewData$downstream_pop_consol[downstream_forNewData$downstream_pop_consol < 0] <- 0
      
      downstream_presence <- ifelse(downstream_forNewData$downstream_pop_consol == 0,0,1)
      newData2[,"downstream_presence"] <- (downstream_presence - mean_downstream)/sd_downstream    
      
      # Recalculate interaction terms
      if (!is.na(interaction_terms)) {
        for (i in 1:length(interaction_terms)) {
          first_set <- interaction_terms[[i]]
          first_term <- newData2[,first_set[1]]
          second_term <- newData2[,first_set[2]]
          newData2[,length(means_relData)+i] <- first_term * second_term
        }
      }
      
      
      # Make new calculations
      attempt <- as.matrix(newData2) %*% as.matrix(beta_mat[s,])
      eta <- sweep(attempt, 2, alpha_mat[s], "+")
      expeta <- exp(eta)
      probabilities <- expeta/(1+expeta)
      
      # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
      new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
      
      # # Introduce upstream lakes again
      # introductions_by_dispersal2 <- introduction_by_dispersal(data_presences$waterBodyID,200, connect_db)
      # new_presences <- ifelse(pop_proximity$new_presences ==1, 1, ifelse(threshold_var > threshold_int, 1, ifelse(pop_proximity$waterBodyID %in% introductions_by_dispersal2,1, 0)))
      
      
      # Establish which lakes now have presences.
      all_presences <- ifelse(pop_proximity$all_presences ==1, 1, ifelse(new_presences == 1, 1, 0))
    }
  }
  time2 <- Sys.time()
  difftime_1 <- round(as.numeric(difftime(time2, time1,
                                          units = "mins")),3)
  if (s %% 20 == 0) {print(paste0("Run ", s, " of ",  n_loops, " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*n_loops/s-difftime_1,5), " minutes left."))}
  periods[,s] <- all_presences
}


print("Finished forecasting.")


forecasts <- periods

if (temp_scenario == TRUE) {
  saveRDS(forecasts, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/forecasts_",population_threshold,"km_wTempIncrease_interaction.RDS"))
} else {
  saveRDS(forecasts, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/forecasts_",population_threshold,"km_interaction.RDS"))
}

