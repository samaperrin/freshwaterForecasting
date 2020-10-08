#---------------------------------------------------------------------#
# 6. Forecasting
#---------------------------------------------------------------------#

focal_species_number <- which(species_list == focal_species)

#---------------------------------------------------------------------#
# 6A. Get Scaling Parameters
#
# Variables within the model will change, but they have to be scaled 
# using the same mands and SDs as our original data, so we need to define
# those first.
#
#---------------------------------------------------------------------#
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

# Import data and draws and define components

model_output <- readRDS(paste0("./Data/whole_model_output_",file_distance_measure,"_pooled.RDS"))
catchment_data <- readRDS("./Data/lakesInCatchments.RDS")
catchment_mapDF <- map_df(catchment_data$catchmentsByLake, ~data.frame(waterBodyID=.),.id="ebint")

# This is the full data. No transformations, and all lakes, even the ones in the native range.
species_preModel_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/species_model_data.RDS"))
species_preModel_data_nonNative <- species_preModel_data %>%
  filter(native == 0)


draws <- model_output$whole_model_params$draws
beta <- model_output$whole_model_params$beta[,focal_species_number]
alpha <- model_output$whole_model_params$alpha[focal_species_number,]

# Whatever data/model we are using for out parameters, we are now applying the parameters to all the data.
# That means we need to scale all the data using whatever means and SDs with which we scaled the data that
# was used for the model. Luckily, we have the means and SDs of our relevant values from a previous script.
standardising_values <- readRDS("./Data/standardising_values.RDS")

means_relData <- standardising_values$means
SDs_relData <- standardising_values$SDs

# This data was collected at line 33-51 of script 4. THis is our environmental data for the model.
env_data <- species_preModel_data_nonNative %>%
  dplyr::transmute(areaL = log(area_km2+1),
                   dist_roadL = log(distance_to_road+1),
                   temp = eurolst_bio10,
                   pop_distL = log(dist_n_pop+1),
                   HFP = HFP,
                   n_pop = log(no_n_pop_20km+1),
                   upstream_presence = upstream_presence,
                   downstream_presence = downstream_presence,
                   nearby_weighted = log(nearby_weightings*10000+1))

env_data_all_almost <- sweep(env_data,2,means_relData,"-")
env_data_all <- sweep(env_data_all_almost,2,SDs_relData,"/")

all_variables <- colnames(env_data)

# Now need to create interaction terms



print("Data is scaled.")

#---------------------------------------------------------------------#
# 6B. Set Figures Beforehand
#
# We need to set up a bunch of parameters to use in the forecasting
# script. 
#
#---------------------------------------------------------------------#

# Get betas and alphas for the looping
beta_init <- as.matrix(calculate(beta,values=draws))
alpha_init <- as.matrix(calculate(alpha,values=draws))

### Get betas and alphas
calc_beta <- apply(beta_init,2,mean)
calc_alpha <- apply(alpha_init,2,mean)

# Get a vector of where we currently have presences outside of the native range
presences <- species_preModel_data_nonNative$introduced


set.seed(147)
chosen_iterations <- sample(1:nrow(alpha_init), n_loops)

beta_mat <- beta_init[chosen_iterations,]
alpha_mat <- alpha_init[chosen_iterations,]

# Need vector of where we have presences outside of the native range at the moment
presence_data <- species_preModel_data[species_preModel_data$presence == 1,]

# These variables define actions in the loop
number_reps <- 5

# These variables will be used to scale projections of numnber of populations nearby and distance to nearest population
mean_weightings <- means_relData["nearby_weighted"]
sd_weightings <- SDs_relData["nearby_weighted"]

mean_upstream <- means_relData["upstream_presence"]
sd_upstream <- SDs_relData["upstream_presence"]

mean_downstream <- means_relData["downstream_presence"]
sd_downstream <- SDs_relData["downstream_presence"]

introductions <- matrix(NA,nrow=nrow(env_data),ncol=1)
periods <- matrix(NA,nrow=nrow(env_data),ncol=n_loops)
# introduction_probs <- list()
# populations <- list()

print("Other stuff in place, starting forecast run.")


time1 <- Sys.time()
plan(multiprocess)
establishment_list_all <- furrr::future_map(.x = 1:n_loops, .progress = TRUE,
                                         .f = ~{
                                           
                                           for (j in 1:number_reps) {
                                             # So, first step is to create predictions based purely off a rise in temperature. Because this is slightly different to what we'll do in the other steps, there needs to be an if clause.
                                             
                                             if (j==1) {
                                               # We create a table containing our new data
                                               
                                               newData <- env_data_all
                                               
                                               # Prepare interaction terms
                                               if (!is.na(interaction_terms)) {
                                                 for (i in 1:length(interaction_terms)) {
                                                   first_set <- interaction_terms[[i]]
                                                   first_term <- newData[,first_set[1]]
                                                   second_term <- newData[,first_set[2]]
                                                   newData[,ncol(env_data_all)+i] <- first_term * second_term
                                                 }
                                               }
                                               
                                               options(warn=-1)
                                               if(!is.na(parameters_to_ignore)) {
                                                 newData <- newData[,-parameters_to_ignore]
                                               } else {newData <- newData}
                                               options(warn=0)
                                               
                                               # Unfortunately because we're predicting for 600k + lakes, it's impossible to run the preditions all at one.
                                               # So the loop below runs them in chunks of 10k lakes at a time.
                                               attempt <- as.matrix(newData) %*% as.matrix(beta_mat[.x,])
                                               eta <- greta::sweep(attempt, 2, alpha_mat[.x], "+")
                                               expeta <- exp(eta)
                                               probabilities <- expeta/(1+expeta)
                                               
                                               # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
                                               new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
                                               all_presences <- ifelse(presences ==1, 1, ifelse(new_presences == 1, 1, 0))
                                               
                                               newData2 <- newData
                                             } else {
                                               
                                               ### And now we move on to the second part of the loop, whereby we redefine those variables pertaining to closest population and numbers of close populations at each step
                                               
                                               # Now figure out the new distance to closest population
                                               pop_proximity <- cbind(all_presences, species_preModel_data_nonNative[,c("utm_x","utm_y","decimalLatitude","decimalLongitude","locationID")])
                                               
                                               # Need to bind our new presences together with all the presences in the native range.
                                               data_presences_nonnative <- as.data.frame(pop_proximity %>% filter(all_presences == 1) %>% 
                                                                                           dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude))
                                               data_presences_native <- species_preModel_data %>% filter(native == 1 & presence == 1) %>% 
                                                 dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude)
                                               data_presences <- rbind(data_presences_nonnative, data_presences_native)
                                               
                                               data_all <- pop_proximity %>%
                                                 dplyr::select(utm_x,utm_y,locationID, decimalLongitude, decimalLatitude) %>%
                                                 distinct() %>%
                                                 as.data.frame()
                                               
                                               
                                               # The get.knnx function returns the distance from a lake in table B to the k closest lakes in table A
                                               
                                               nn <- get.knnx(data_presences[c("utm_x","utm_y")],data_all[c("utm_x","utm_y")],k=100)
                                               
                                               weighting_nearby_pop <- nn$nn.dist
                                               weighting_nearby_pop[weighting_nearby_pop > 100000 | weighting_nearby_pop == 0] <- NA
                                               idw <- 1/weighting_nearby_pop
                                               nearby_inv_weightings_vec <- rowSums(idw,na.rm=TRUE)
                                               
                                               nearby_inv_weightings_table <- cbind(nearby_inv_weightings_vec,data_all$locationID)
                                               colnames(nearby_inv_weightings_table) <- c("weightings", "locationID")
                                               
                                               ordered_weightings <- merge(species_preModel_data_nonNative["locationID"],nearby_inv_weightings_table,all.x=TRUE,by="locationID")
                                               matched_ordered_weightings <- ordered_weightings[match(species_preModel_data_nonNative[,"locationID"],ordered_weightings$locationID),]  
                                               
                                               
                                               new_weightings <- log(as.numeric(as.character(matched_ordered_weightings$weightings))*10000+1)
                                               
                                               # So now that we have the new measurements for closest population, these need to be scaled against those that we had for the initial population. So we subtract the mean and divide by the standard deviation.
                                               newData2[,"nearby_weighted"] <- (new_weightings - mean_weightings)/sd_weightings
                                               
                                               # Now we need the number of populations upstream
                                               catchment_mapDF_ts <- catchment_mapDF
                                               
                                               # Get list of all water bodies with presence
                                               all_presence_table <- species_preModel_data[,c("locationID","waterBodyID")]
                                               all_presence_table$present <- ifelse(all_presence_table$locationID %in% data_presences$locationID, 1, 0)
                                               
                                               all_WB_present <- all_presence_table[all_presence_table$present == 1,"waterBodyID"]
                                               catchment_mapDF_ts$present <- ifelse(catchment_mapDF_ts$waterBodyID %in% all_WB_present, 1,0)
                                               # Define how many populations there are in each catchment
                                               new_upstream_pops <- catchment_mapDF_ts %>%
                                                 group_by(ebint) %>%
                                                 summarize(upstream_pops = sum(present))
                                               upstream_forNewData <- merge(species_preModel_data[species_preModel_data$native==0,
                                                                                                  c("waterBodyID","ebint","locationID")], 
                                                                            new_upstream_pops, all.x=TRUE, by = "ebint")
                                               matched_upstream_forNewData <- upstream_forNewData[match(species_preModel_data_nonNative[,"locationID"],upstream_forNewData$locationID),]
                                               # We still have NAs where a lake isn't in its corresponding catchment
                                               # These can become zeroes, as there clearly isn't anything upstream in that catchment anyway
                                               matched_upstream_forNewData$upstream_pops[is.na(matched_upstream_forNewData$upstream_pops)] <- 0
                                               # If we have a situation where we have 2 presences and one of them is the focal lake, we need to reduce
                                               # the number of upstream lakes to 1.
                                               matched_upstream_forNewData$upstream_pop_consol <- ifelse(matched_upstream_forNewData$waterBodyID %in% all_WB_present,
                                                                                                         matched_upstream_forNewData$upstream_pops - 1, matched_upstream_forNewData$upstream_pops)
                                               # We may have a couple of stray negatives. These can simply become 0.
                                               matched_upstream_forNewData$upstream_pop_consol[matched_upstream_forNewData$upstream_pop_consol<0] <- 0
                                               
                                               upstream_presence <- ifelse(matched_upstream_forNewData$upstream_pop_consol == 0,0,1)
                                               newData2[,"upstream_presence"] <- (upstream_presence - mean_upstream)/sd_upstream
                                               
                                               # Need the number of populations downstream now
                                               basin_mapDF <- species_preModel_data[,c("waterBodyID","eb_waterregionID")] 
                                               basin_mapDF <- basin_mapDF[!duplicated(basin_mapDF),]
                                               basin_mapDF$present <- ifelse(basin_mapDF$waterBodyID %in% all_WB_present, 1,0)
                                               # Define how many populations there are in each basin
                                               new_downstream_pops <- basin_mapDF %>%
                                                 group_by(eb_waterregionID) %>%
                                                 summarize(downstream_pops = sum(present))
                                               downstream_forNewData <- merge(species_preModel_data_nonNative[,c("waterBodyID","eb_waterregionID","locationID")], new_downstream_pops
                                                                              , all.x=TRUE, by = "eb_waterregionID")
                                               
                                               matched_downstream_forNewData <- downstream_forNewData[match(species_preModel_data_nonNative[,"locationID"],downstream_forNewData$locationID),]
                                               matched_downstream_forNewData$downstream_pop_consol <- ifelse(matched_downstream_forNewData$waterBodyID %in% all_WB_present,
                                                                                                             matched_downstream_forNewData$downstream_pops - 1, matched_downstream_forNewData$downstream_pops)
                                               matched_downstream_forNewData$downstream_pop_consol <- matched_downstream_forNewData$downstream_pop_consol - matched_upstream_forNewData$upstream_pop_consol
                                               # Once more, get rid of negatives
                                               matched_downstream_forNewData$downstream_pop_consol[matched_downstream_forNewData$downstream_pop_consol < 0] <- 0
                                               
                                               downstream_presence <- ifelse(matched_downstream_forNewData$downstream_pop_consol == 0,0,1)
                                               newData2[,"downstream_presence"] <- (downstream_presence - mean_downstream)/sd_downstream    
                                               
                                               # Recalculate interaction term
                                               if (!is.na(interaction_terms)) {
                                                 for (i in 1:length(interaction_terms)) {
                                                   first_set <- interaction_terms[[i]]
                                                   first_term <- newData2[,all_variables[first_set[1]]]
                                                   second_term <- newData2[,all_variables[first_set[2]]]
                                                   newData2[,ncol(env_data_all)-length(parameters_to_ignore)+i] <- first_term * second_term
                                                 }
                                               }
                                               
                                               # Make new calculations
                                               attempt <- as.matrix(newData2) %*% as.matrix(beta_mat[.x,])
                                               eta <- sweep(attempt, 2, alpha_mat[.x], "+")
                                               expeta <- exp(eta)
                                               probabilities <- expeta/(1+expeta)
                                               
                                               # Establish which lakes now have presences by asking whether or not we have an introduction, based on a bernoulli estimate.
                                               new_presences <- rbinom(length(probabilities), size = 1, prob=probabilities)
                                               all_presences <- ifelse(all_presences ==1, 1, ifelse(new_presences == 1, 1, 0))
                                             }
                                           }
                                           all_presences
                                         })
time2 <- Sys.time()
print(time2-time1)

introductions <- matrix(NA, nrow=length(establishment_list_all[[1]]), ncol = n_loops)
for (l in 1:n_loops) {
  introductions[,l] <- establishment_list_all[[l]]
}


saveRDS(introductions, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/forecasts_concurrent_",file_distance_measure,".RDS"))




