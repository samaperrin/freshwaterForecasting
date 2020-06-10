#---------------------------------------------------------------------#
# 6. Historical Forecasting
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 6A. Get Historical Data
#
# We need to organsie the data so it looks the way it did in 1969. That
# means turning introductions after that eyar into absences, as well
# as redoing all parameters which rely on number of upstream and 
# downstream populations.
#
#---------------------------------------------------------------------#

model_output <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                               "/whole_model_output_",population_threshold,"km.RDS"))
model_data_extra <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                   "/whole_model_data_",population_threshold,"km.RDS"))
species_model_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                     "/species_model_data.RDS"))
species_year_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                    "/species_byYear.RDS"))

catchment_data <- readRDS("./Data/lakesInCatchments.RDS")
catchment_mapDF <- map_df(catchment_data$catchmentsByLake, ~data.frame(waterBodyID=.),.id="ebint")

# Merge model data with year data
model_data_withYears <- merge(model_data_extra$raw_data, species_year_data, all.x=TRUE, by.x="no_vatn_lnr",by.y = "vatnLnr")

# Get rid of all introductions without a year present
model_data_withYears <- model_data_withYears %>% filter(!(is.na(firstYear) & introduced == 1))

# Create our paramater for whether or not the species was present back then
model_data_withYears$historicIntroductions <- ifelse(is.na(model_data_withYears$firstYear),0,
                                                     ifelse(model_data_withYears$firstYear < 1969, 1, 0))
model_data_withYears$historicPresences <- ifelse(model_data_withYears$presence == 0, 0, 
                                                 ifelse(model_data_withYears$native == 1, 1, 
                                                        ifelse(model_data_withYears$firstYear > 1969, 0,1)))
  
# Now make every leftover presence's year 0
model_data_withYears$firstYear <- ifelse(is.na(model_data_withYears$firstYear), 0 ,model_data_withYears$firstYear)

#### Now we create data for the 1969 stuff #####

# Firstly we produce a data frame of our introductions.
species_intro_data <- model_data_withYears %>% 
  filter(historicIntroductions == 1)

# We also produce a dataset with all presences. These won't be in the model, but they'll be used to calculate nearest population
# and number of nearby populations.
species_prese_data <- model_data_withYears %>%
  filter(historicPresences == 1)

# We also need to produce a series of absences.
species_absen_data <- model_data_withYears %>%
  filter(native == 0 & historicPresences == 0)

# We calculate nearby populations for the absences first. We take all current presences into account here,
# there's no need to put in a specific year.
# Calculate nearest population first.
nn_nearest_abs <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=1)
species_absen_data$dist_n_pop <- nn_nearest_abs$nn.dist[,1]

# THen get all nearby populations.
k_abs <- ifelse(nrow(species_prese_data)<70,nrow(species_prese_data),70)
nn_nearby_abs <- get.knnx(species_prese_data[c("utm_x","utm_y")],species_absen_data[c("utm_x","utm_y")],k=k_abs)

# File this down to all populations within chosen radius.
number_nearby_pop <- nn_nearby_abs$nn.dist
number_nearby_pop[number_nearby_pop < population_threshold*1000] <- 1
number_nearby_pop[number_nearby_pop >= population_threshold*1000] <- 0
number_nearby_pop <- rowSums(number_nearby_pop)

# We can now do this for populations upstream
# Need to group all lakes by the catchments they're in
catchmentsByLake <- catchment_data$catchmentsByLake
catchment_mapDF <- map_df(catchmentsByLake, ~data.frame(waterBodyID=.),.id="ebint")
catchment_mapDF_abs <- catchment_mapDF

# Turn lakes with presences in catchment to 1s, absences to 0s
catchment_mapDF_abs$waterBodyID[!(catchment_mapDF_abs$waterBodyID %in% species_prese_data$waterBodyID)] <- 0
catchment_mapDF_abs$waterBodyID[catchment_mapDF_abs$waterBodyID %in% species_prese_data$waterBodyID] <- 1

# And get number of populations
upstream_pops <- catchment_mapDF_abs %>%
  group_by(ebint) %>%
  summarize(upstream_pops = sum(waterBodyID))

# And for populations downstream we do exactly the same thing
basin_mapDF <- model_data_withYears[,c("waterBodyID","eb_waterregionID")]
basin_mapDF_abs <- basin_mapDF

basin_mapDF_abs$waterBodyID[!(basin_mapDF_abs$waterBodyID %in% species_prese_data$waterBodyID)] <- 0
basin_mapDF_abs$waterBodyID[basin_mapDF_abs$waterBodyID %in% species_prese_data$waterBodyID] <- 1

basin_pops <- basin_mapDF_abs %>%
  group_by(eb_waterregionID) %>%
  summarize(all_pops = sum(waterBodyID))

upstream_pops_merged <- merge(species_absen_data[,c("ebint","waterBodyID","eb_waterregionID")],upstream_pops,all.x=TRUE,by="ebint")
all_pops_merged <- merge(upstream_pops_merged, basin_pops,all.x=TRUE, by="eb_waterregionID")

# You might notice at this point that we have NAs in our upstream pops section. THis is because
# some lakes fall outside of the geometries constructed for their catchments and the 
# catchment therefore has no lake inside it. This only occurs for catchments with one lake,
# which means they're the futhest upstream lake anyway, so we can just turn the NAs to 0s.
all_pops_merged[is.na(all_pops_merged$upstream_pops),] <- 0

all_pops_merged$downstream_pops <- all_pops_merged$all_pops - all_pops_merged$upstream_pops

# Correct any mess ups to zero
all_pops_merged[all_pops_merged$downstream_pops<0,] <- 0


# Put this into a data frame which we will then add more to later on.
spatial_data <- data.frame(species_absen_data$no_vatn_lnr, nn_nearest_abs$nn.dist, 
                           number_nearby_pop, all_pops_merged$upstream_pops, all_pops_merged$downstream_pops)
colnames(spatial_data) <- c("no_vatn_lnr","dist_n_pop","no_n_pop",
                            "upstream_pops","downstream_pops")

# Now we cycle through year by year. All presences in native range are assumed to have been there ebfore the
# first time-step.
time_steps <- sort(unique(species_intro_data$firstYear)) 
for (i in 1:length(time_steps)) {
  year_step <- time_steps[i]
  species_intro_historic <- species_intro_data %>% filter(firstYear == year_step)
  species_presence_historic <- species_prese_data %>% filter(firstYear < year_step)
  
  # Special case for first time step for Rainbow Trout, since they have no historic range in Norway.
  if (nrow(species_presence_historic) == 0) {
    spatial_data_timeStep <- data.frame(species_intro_historic$no_vatn_lnr, 9999999, 0, 0, 0)
    colnames(spatial_data_timeStep) <- c("no_vatn_lnr","dist_n_pop","no_n_pop",
                                         "upstream_pops","downstream_pops")
  } else {
    
    nn_nearest <- get.knnx(species_presence_historic[c("utm_x","utm_y")],species_intro_historic[c("utm_x","utm_y")],k=2)
    # Same as above, but need to make sure we don't include a lake as its own enarest population
    nearest_pop_vec <- ifelse(nn_nearest$nn.dist[,1]==0,nn_nearest$nn.dist[,2],nn_nearest$nn.dist[,1])
    
    k_nearby <- ifelse(nrow(species_presence_historic) < 100, nrow(species_presence_historic)-1, 100)
    
    nn_nearby <- get.knnx(species_presence_historic[c("utm_x","utm_y")],species_intro_historic[c("utm_x","utm_y")],k=k_nearby)
    
    # File this down to all populations within chosen radius.
    number_nearby_pop <- nn_nearby$nn.dist
    number_nearby_pop[number_nearby_pop < population_threshold*1000] <- 1
    number_nearby_pop[number_nearby_pop >= population_threshold*1000] <- 0
    number_nearby_pop <- rowSums(number_nearby_pop)
    
    # Get populations upstream
    catchment_mapDF_tS <- catchment_mapDF
    
    catchment_mapDF_tS$waterBodyID[!(catchment_mapDF_tS$waterBodyID %in% species_presence_historic$waterBodyID)] <- 0
    catchment_mapDF_tS$waterBodyID[catchment_mapDF_tS$waterBodyID %in% species_presence_historic$waterBodyID] <- 1
    
    upstream_pops <- catchment_mapDF_tS %>%
      group_by(ebint) %>%
      summarize(upstream_pops = sum(waterBodyID))
    
    # And populations downstream
    basin_mapDF_tS <- basin_mapDF
    
    basin_mapDF_tS$waterBodyID[!(basin_mapDF_tS$waterBodyID %in% species_presence_historic$waterBodyID)] <- 0
    basin_mapDF_tS$waterBodyID[basin_mapDF_tS$waterBodyID %in% species_presence_historic$waterBodyID] <- 1
    
    basin_pops <- basin_mapDF_abs %>%
      group_by(eb_waterregionID) %>%
      summarize(all_pops = sum(waterBodyID))
    
    # Bring them together so we can get number of downstream populations
    upstream_pops_merged <- merge(species_intro_historic[,c("ebint","waterBodyID","eb_waterregionID")],upstream_pops,all.x=TRUE,by="ebint")
    all_pops_merged <- merge(upstream_pops_merged, basin_pops,all.x=TRUE, by="eb_waterregionID")
    all_pops_merged[is.na(all_pops_merged$upstream_pops),] <- 0
    
    all_pops_merged$downstream_pops <- all_pops_merged$all_pops - all_pops_merged$upstream_pops
    all_pops_merged[all_pops_merged$downstream_pops<0,] <- 0
    
    
    spatial_data_timeStep <- data.frame(species_intro_historic$no_vatn_lnr, nearest_pop_vec, 
                                        number_nearby_pop, all_pops_merged$upstream_pops, all_pops_merged$downstream_pops)
    colnames(spatial_data_timeStep) <- c("no_vatn_lnr","dist_n_pop","no_n_pop",
                                         "upstream_pops","downstream_pops")
  }
  
  spatial_data <- rbind(spatial_data,spatial_data_timeStep)
  
}

# Now, get upstream/downstream populations
spatial_data$upstream_presence <- ifelse(spatial_data$upstream_pops == 0, 0, 1)
spatial_data$downstream_presence <- ifelse(spatial_data$downstream_pops == 0, 0, 1)

#---------------------------------------------------------------------#
# 6B. Prepare data for forecasting
#
# Need to get everything ready to run through the forecasting loop.

#
#---------------------------------------------------------------------#


# Now that we have these variables, need to produce a data frame free of incomplete cases.
# HOWEVER, we still need lakes outside the native range to find total number of rpesences
# nearby, so we keep these
base_data <- model_data_withYears[,c(1:20,30:32)]


historic_all_data <- merge(base_data,spatial_data, all.x=TRUE, by="no_vatn_lnr")

# Get the model stuff in too
draws <- model_output$draws
beta <- model_output$beta
p <- model_output$p
alpha <- model_output$alpha


# So now we have data as it looked in 1969. We need to scale it via the model parameters.
model_data <- model_data_extra$raw_data %>%
  filter(native == 0)

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

historic_model_data <- historic_all_data %>%
  filter(native==0)

historic_envData <- historic_model_data %>%
  transmute(areaL = log(area_km2+1),
            dist_roadL = log(distance_to_road+1),
            temperature = eurolst_bio10,
            pop_distL = log(dist_n_pop+1),
            HFP = HFP,
            n_pop = log(no_n_pop+1),
            upstream_presence = upstream_presence,
            downstream_presence = downstream_presence)

scaled_envData_withYears1 <- greta::sweep(historic_envData, 2, means_relData, "-")
historic_envData <- greta::sweep(scaled_envData_withYears1, 2, SDs_relData, "/")

# Add in our two interaction terms
if (!is.na(interaction_terms)) {
  for (i in 1:length(interaction_terms)) {
    first_set <- interaction_terms[[i]]
    first_term <- historic_envData[,first_set[1]]
    second_term <- historic_envData[,first_set[2]]
    if (i == 1) {interaction_term <- first_term * second_term
    } else {
      interaction_term <- cbind(interaction_term, first_term * second_term)
    }
  }
  historic_envData <- cbind(historic_envData, interaction_term)
}

### Get betas and alphas
calc_beta <- apply(as.matrix(calculate(beta,values=draws)),2,mean)
calc_alpha <- apply(as.matrix(calculate(alpha,values=draws)),2,mean)

# Get betas and alphas for the looping
beta_init <- as.matrix(calculate(beta,values=draws))
alpha_init <- as.matrix(calculate(alpha,values=draws))

set.seed(147)
n_loops_historic <- n_loops
chosen_iterations <- sample(1:nrow(alpha_init), n_loops_historic)

beta_mat <- beta_init[chosen_iterations,]
alpha_mat <- alpha_init[chosen_iterations,]

# Need vector of where we have presences outside of the native range at the moment
presences <- historic_model_data$historicIntroductions

nn_all <- get.knnx(historic_all_data[,c("utm_x","utm_y")],historic_model_data[c("utm_x","utm_y")],k=100)
nn_all$nn.index[nn_all$nn.dist > population_threshold * 1000 | nn_all$nn.dist == 0] <- 0


# These variables define actions in the loop
number_reps <- 5

# This defines in incremental increase over 50 years based on a 1 degree increase in temperature over 50 years.
# Still need to figure out what to use here!
temp_increase <- 0 # This needs to be replaced with whatever the increase in average temperature was between
# the two timepoints (1969-2019)
temp_step <- ifelse(temp_scenario == TRUE,temp_increase/SDs_relData["temp"],rep(0,length(SDs_relData["temp"])))

# These variables will be used to scale projections of numnber of populations nearby and distance to nearest population
mean_n_pop <- means_relData["n_pop"]
sd_n_pop <- SDs_relData["n_pop"]

sd_dist <- SDs_relData["pop_distL"]
mean_dist <- means_relData["pop_distL"]

mean_upstream <- means_relData["upstream_presence"]
sd_upstream <- SDs_relData["upstream_presence"]

mean_downstream <- means_relData["downstream_presence"]
sd_downstream <- SDs_relData["downstream_presence"]

introductions <- matrix(NA,nrow=nrow(historic_envData),ncol=1)
periods <- matrix(NA,nrow=nrow(historic_envData),ncol=n_loops_historic)
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
for (s in 1:n_loops_historic) {
  for (j in 1:number_reps) {
    # So, first step is to create predictions based purely off a rise in temperature. Because this is slightly different to what we'll do in the other steps, there needs to be an if clause.
    
    if (j==1) {
      # We create a table containing our new data
      
      tempIncrease <- rnorm(nrow(historic_envData), temp_step/5, 0.05)
      newData <- historic_envData
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
      pop_proximity <- cbind(all_presences, historic_model_data[,c("utm_x","utm_y","decimalLatitude","decimalLongitude","locationID")])
      
      # Need to bind our new presences together with all the presences in the native range.
      data_presences_nonnative <- as.data.frame(pop_proximity %>% filter(all_presences == 1) %>% 
                                                  dplyr::select(locationID, utm_x, utm_y, decimalLongitude, decimalLatitude))
      data_presences_native <- historic_all_data %>% filter(native == 1 & presence == 1) %>% 
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
      
      ordered_distances <- merge(historic_model_data,distance_data,all.x=TRUE,by="locationID") %>%
        arrange(no_vatn_lnr)
      
      new_distances <- log(as.numeric(as.character(ordered_distances$dist_to_closest_pop))+1)
      
      # So now that we have the new measurements for closest population, these need to be scaled against those that we had for the initial population. So we subtract the mean and divide by the standard deviation.
      newData2[,"pop_distL"] <- (new_distances-mean_dist)/sd_dist
      
      # Now we need to get the new measurements for number of close populations
      all_presence_table <- historic_all_data[,c("locationID","waterBodyID","no_vatn_lnr")]
      all_presence_table$present <- ifelse(all_presence_table$locationID %in% data_presences$locationID, 1, 0)
      
      rowNumbers_withPresence <- which(all_presence_table$present==1)
      nn_inds <- nn_all$nn.index
      nn_inds[!(nn_inds %in% rowNumbers_withPresence)] <- 0
      nn_inds[nn_inds %in% rowNumbers_withPresence] <- 1
      number_nearby_pop_vec <- log(rowSums(nn_inds)+1)
      # 
      # Scale it like we did for the last variable
      
      newData2[,"n_pop"] <- (number_nearby_pop_vec - mean_n_pop)/sd_n_pop
      
      # Now we need the number of populations upstream
      catchment_mapDF_ts <- catchment_mapDF
      
      # Get list of all water bodies with presence
      all_WB_present <- all_presence_table[all_presence_table$present == 1,"waterBodyID"]
      catchment_mapDF_ts$present <- ifelse(catchment_mapDF_ts$waterBodyID %in% all_WB_present, 1,0)
      # Define how many populations there are in each catchment
      new_upstream_pops <- catchment_mapDF_ts %>%
        group_by(ebint) %>%
        summarize(upstream_pops = sum(present))
      upstream_forNewData <- merge(historic_model_data[,c("waterBodyID","ebint","no_vatn_lnr")], new_upstream_pops, all.x=TRUE, by = "ebint") %>%
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
      basin_mapDF <- historic_all_data[,c("waterBodyID","eb_waterregionID")] 
      basin_mapDF <- basin_mapDF[!duplicated(basin_mapDF),]
      basin_mapDF$present <- ifelse(basin_mapDF$waterBodyID %in% all_WB_present, 1,0)
      # Define how many populations there are in each basin
      new_downstream_pops <- basin_mapDF %>%
        group_by(eb_waterregionID) %>%
        summarize(downstream_pops = sum(present))
      downstream_forNewData <- merge(historic_model_data[,c("waterBodyID","eb_waterregionID","no_vatn_lnr")], new_downstream_pops
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
  if (s %% 20 == 0) {print(paste0("Run ", s, " of ",  n_loops_historic, " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*n_loops_historic/s-difftime_1,5), " minutes left."))}
  periods[,s] <- all_presences
}


print("Finished forecasting. Downloading map of Norway.")

#---------------------------------------------------------------------#
# 6D. Simple visualisation
#
# Compare forecast to what actually happened.

#
#---------------------------------------------------------------------#

lake_likelihoods_noIncrease <- apply(periods, 1, mean)

all_data_likelihoods <- cbind(historic_model_data,lake_likelihoods_noIncrease)


# Produce a map of Norway to use
Norway<-getData("GADM", country="NO", level=0)
Norway1<-getData("GADM", country="NO", level=1)

# Need to get it down to just the southern counties
Norway1_sub<-Norway1[!(Norway1@data$NAME_1 %in% c('Troms', 'Finnmark', 'Nordland')),]
par(mar=c(1,1,1,1))
Norway_df <- fortify(Norway)
Norway1_sub_df <- fortify(Norway1_sub)
Norway1_sub_sf <- st_as_sf(Norway1_sub)

Norway_asOne <- st_union(Norway1_sub_sf)

print("Downloading species native range")

# Get native distribution for everythign except Rainbow Trout
if(focal_species != "Oncorhynchus mykiss") {
  hk_distribution_map <- readRDS("Data/native_distribution.rds")
  species_native <- hk_distribution_map[hk_distribution_map$canonicaln == focal_species & 
                                          hk_distribution_map$establishm == "native",]
  species_native_intersection <- st_intersection(st_make_valid(species_native), st_as_sf(Norway1_sub))
  species_native_intersection_union <- st_union(species_native_intersection)
}

predicted_appearances_historic <- ggplot(data=all_data_likelihoods) + 
  geom_sf(data = Norway_asOne, fill = 'white',linetype=3) +
  stat_summary_hex(aes(x = decimalLongitude, 
                       y = decimalLatitude,
                       z = lake_likelihoods_noIncrease),
                   bins=120) +
  scale_fill_gradient2(low = 'light blue', mid="yellow", high="red", midpoint = 0.5,
                       breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1)) +
  #geom_hex(data = filter(all_data_likelihoods, introduced == 1) ,
  #         aes(x = decimalLongitude, 
  #             y = decimalLatitude),
  #        fill="green",
  #         bins=120) +
  geom_hex(data = filter(all_data_likelihoods, historicIntroductions == 1) ,
           aes(x = decimalLongitude, 
               y = decimalLatitude),
           fill="black",
           bins=120) +
  labs(subtitle=paste0("Forecasted distribution of ",focal_species," in 50 Years"), fill='Establishment risk') +
  xlab("Longitude") +
  ylab("Latitude")

if (focal_species != "Oncorhynchus mykiss") {
  predicted_appearances_historic <- predicted_appearances_historic + 
    geom_sf(data = species_native_intersection_union, fill = "dark grey")
}


predicted_appearances_historic_wPresent <- ggplot(data=all_data_likelihoods) + 
  geom_sf(data = Norway_asOne, fill = 'white',linetype=3) +
  stat_summary_hex(aes(x = decimalLongitude, 
                       y = decimalLatitude,
                       z = lake_likelihoods_noIncrease),
                   bins=120) +
  scale_fill_gradient2(low = 'light blue', mid="yellow", high="red", midpoint = 0.5,
                       breaks=seq(0,1,0.25), #breaks in the scale bar
                       limits=c(0,1)) +
  geom_hex(data = filter(all_data_likelihoods, introduced == 1) ,
          aes(x = decimalLongitude,
              y = decimalLatitude),
         fill="green",
          bins=120) +
  geom_hex(data = filter(all_data_likelihoods, historicIntroductions == 1) ,
           aes(x = decimalLongitude, 
               y = decimalLatitude),
           fill="black",
           bins=120) +
  labs(subtitle=paste0("Forecasted distribution of ",focal_species," in 50 Years"), fill='Establishment risk') +
  xlab("Longitude") +
  ylab("Latitude")

if (focal_species != "Oncorhynchus mykiss") {
  predicted_appearances_historic_wPresent <- predicted_appearances_historic_wPresent + 
    geom_sf(data = species_native_intersection_union, fill = "dark grey")
}

maps <- list(map_onlyHistoric = predicted_appearances_historic, map_withPresent = predicted_appearances_historic_wPresent)

all_data_likelihoods$new_introductions <- ifelse(all_data_likelihoods$introduced == 1 & all_data_likelihoods$historicIntroductions == 0,
                                                 1,0)

historic_forecast <- list(historic_forecast = all_data_likelihoods, historic_maps = maps)
saveRDS(historic_forecast, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/forecasts_",population_threshold,"km_historic.RDS"))

summary(all_data_likelihoods[all_data_likelihoods$new_introductions==1,"lake_likelihoods_noIncrease"])
summary(all_data_likelihoods[,"lake_likelihoods_noIncrease"])






