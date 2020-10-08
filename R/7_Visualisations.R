#---------------------------------------------------------------------#
# 7. Visualise results of the forecasting
#---------------------------------------------------------------------#

# This script visualises where species are likely to be found after a 2.1 degree increase in
# temperature. 

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

forecasts_fullAbs <- readRDS(file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                            "/forecasts_concurrent_weightedDist_fullAbs.RDS"))

forecasts_pseudoAbs <- readRDS(file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),
                                         "/forecasts_concurrent_weightedDist_pseudoAbs.RDS"))

#---------------------------------------------------------------------#
# 7A. Create table of forecasting stats
# 
# Basically we want to create an enormous table with introduction
# probabilities, uncertainties, effects of temperature.
#------------------------------------------------------------------ ---#

lake_likelihoods_fullAbs <- apply(forecasts_fullAbs, 1, mean)
lake_likelihoods_pseudoAbs <- apply(forecasts_pseudoAbs, 1, mean)

all_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_data_weightedDist.RDS"))
analytics <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/whole_model_analytics_weightedDistinteraction.RDS"))

model_data <- all_data$raw_data %>% filter(native == 0)

# # Get uncertainty variables
# Uncertainty <- analytics$intervals$width
# InitProbs <- analytics$intervals$mean

# Put together one data frame to contain all predictive information
all_data_likelihoods <- cbind(model_data,lake_likelihoods_fullAbs, lake_likelihoods_pseudoAbs)


# Save data for use in QGIS
saveRDS(all_data_likelihoods, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/introduction_likelihoods_weightedDist.RDS"))
write.csv(all_data_likelihoods, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/introduction_likelihoods_weightedDist.csv"))

print("Saved mapping data, moving on to introductions maps")

#---------------------------------------------------------------------#
# 7B. Start mapping
#
# The first of these maps will show average likelihood of establishment
# in each hexagon for each species. Hexagons in which the species was
# already established before the scenario ran are marked black.
# Species native distributions are denoted by the grey polygons.
#
#---------------------------------------------------------------------#

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

# Get native distribution for everythign except Rainbow Trout

  hk_distribution_map <- readRDS("Data/native_distribution.rds")
  species_native <- hk_distribution_map[hk_distribution_map$canonicaln == focal_species & 
                                          hk_distribution_map$establishm == "native",]
  species_native_intersection <- st_intersection(st_buffer(species_native,0), st_as_sf(Norway1_sub))
  species_native_intersection_union <- st_union(species_native_intersection)



# Let's first plot the initial appearances

  # Now let's plot predicted appearances without temp increase, using likelihood of introduction 
  # as our value to average over hexagons
  predicted_appearances_fullAbs <- ggplot(data=all_data_likelihoods) + 
    geom_sf(data = Norway_asOne, fill = 'white',linetype=3) +
    stat_summary_hex(aes(x = decimalLongitude, 
                         y = decimalLatitude,
                         z = lake_likelihoods_fullAbs),
                     bins=n_bins) +
    scale_fill_gradient2(low = 'light blue', mid="yellow", high="red", midpoint = 0.5,
                         breaks=seq(0,1,0.25), #breaks in the scale bar
                         limits=c(0,1)) +
    geom_hex(data = filter(all_data_likelihoods, introduced == 1) ,
             aes(x = decimalLongitude, 
                 y = decimalLatitude),
             fill="black",
             bins=n_bins) + 
    geom_sf(data = species_native_intersection_union, fill = "dark grey") +
    labs(subtitle=paste0(common_names[which(species_list == focal_species)]), fill='Establishment risk') +
    xlab("Longitude") +
    ylab("Latitude")
  
  predicted_appearances_pseudoAbs <- ggplot(data=all_data_likelihoods) + 
    geom_sf(data = Norway_asOne, fill = 'white',linetype=3) +
    stat_summary_hex(aes(x = decimalLongitude, 
                         y = decimalLatitude,
                         z = lake_likelihoods_pseudoAbs),
                     bins=n_bins) +
    scale_fill_gradient2(low = 'light blue', mid="yellow", high="red", midpoint = 0.5,
                         breaks=seq(0,1,0.25), #breaks in the scale bar
                         limits=c(0,1)) +
    geom_hex(data = filter(all_data_likelihoods, introduced == 1) ,
             aes(x = decimalLongitude, 
                 y = decimalLatitude),
             fill="black",
             bins=n_bins) + 
    geom_sf(data = species_native_intersection_union, fill = "dark grey") +
    labs(subtitle=paste0(common_names[which(species_list == focal_species)]), fill='Establishment risk') +
    xlab("Longitude") +
    ylab("Latitude")
  
  if (scale_bar == TRUE) {
    predicted_appearances_fullAbs <- predicted_appearances_fullAbs + 
      ggsn::scalebar(dist = 100, dist_unit = "km",
                     transform = TRUE, model = "WGS84", height =0.02,
                     st.size = 3, st.dist = 0.02, anchor = c(x = 13,y = 58), 
                     x.min = 12, x.max = 14, y.min = 57.99, y.max = 65.12616, 
                     st.bottom=TRUE, border.size = 0.5)
    predicted_appearances_pseudoAbs <- predicted_appearances_pseudoAbs + 
      ggsn::scalebar(dist = 100, dist_unit = "km",
                     transform = TRUE, model = "WGS84", height =0.02,
                     st.size = 3, st.dist = 0.02, anchor = c(x = 13,y = 58), 
                     x.min = 12, x.max = 14, y.min = 57.99, y.max = 65.12616, 
                     st.bottom=TRUE, border.size = 0.5)
  }
  
  if (compass == TRUE) {
    predicted_appearances_fullAbs <- north2(predicted_appearances_fullAbs, x=0.5,y=0.8, scale = 0.1, symbol=2)
    predicted_appearances_pseudoAbs <- north2(predicted_appearances_pseudoAbs, x=0.2,y=0.8, scale = 0.1, symbol=2)  
  }
  
#---------------------------------------------------------------------#
# 7C. Map water basin averages
#
# These maps will now show the average likelihood of establishment 
# across entire water basins, excluding lakes which already had 
# establishments before the scenarios were run. As with the above maps,
# the second graph shows changes with introduction of a temperature
# change.
#
#---------------------------------------------------------------------#

  if (give_basin_plots == TRUE) {
    
    # Now let's plot averages for water basins
    if (download_basins == TRUE) {
      pg_user=rstudioapi::askForPassword("Wallace username")
      pg_password=rstudioapi::askForPassword("password")
      
      # Need two more connections to get the catchment and basin geometries
      pool_basin <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
                           host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
                           idleTimeout = 36000000,
                           options="-c search_path=\"Hydrography\""
      )
      con_basin <- poolCheckout(pool_basin)
      
      # Download basin geometries
      basin_import <- tbl(con_basin,"waterregions_dem_10m_nosefi") %>%
        dplyr::select(gid, geom) %>%
        filter(gid %in% !!all_data_likelihoods$eb_waterregionID) %>%
        collect()
      basin_list <- vector()
      
      # Convert text strings into geometries
      for(i in seq_len(nrow(basin_import))) {
        wkb_int <- structure(list(basin_import$geom[i]), class = "WKB")
        basin_list <- c(basin_list,wkb_int)
      }
      basin_geom <- st_as_sfc(basin_list, EWKB = TRUE)
      basin_geom_sf <- basin_geom %>%
        st_sf %>%
        st_cast
      
      # Give each object a name and save to file
      basin_geom_sf$eb_waterregionID <- basin_import$gid
      saveRDS(basin_geom_sf,file="./Data/waterRegions.rds")
    } else {
      basin_geom_sf <- readRDS("./Data/waterRegions.rds")
    }
    
    # Get just the lakes which didn't already have presences
    all_data_likelihoods_basinTransform <- all_data_likelihoods %>%
      filter(presence == 0)
    
    # Find average likelihoods for basins
    basin_aggregate <- all_data_likelihoods %>%
      group_by(eb_waterregionID) %>%
      summarise(avg_pseudoAbs = mean(lake_likelihoods_pseudoAbs),
                avg_fullAbs = mean(lake_likelihoods_fullAbs))
    
    # Match up the crs's
    basin_geom_sf <- st_transform(basin_geom_sf, 4326)
    
    basin_plot_data <- merge(basin_geom_sf,basin_aggregate,all.x=TRUE,by="eb_waterregionID")
    # basin_plot_data$average <- ifelse(basin_plot_data$present == 1 | basin_plot_data$native == 1, 1, basin_plot_data$avg_noIncrease)
    
    # And plot
    basin_map_noIncrease <- ggplot(data=basin_plot_data) +
      geom_sf(mapping = aes(fill=avg_pseudoAbs)) +
      scale_fill_gradient2(low = 'light blue', mid="yellow", high="red", midpoint=0.5) +
      labs(subtitle=paste0("Likelihood of introduction of ",focal_species," in 50 Years by drainage basin"), fill='Establishment risk') +
      xlab("Longitude") +
      ylab("Latitude")
    
    # And plot with temperature increase
    basin_map_wIncrease <- ggplot(data=basin_plot_data) +
      geom_sf(mapping = aes(fill=avg_fullAbs)) +
      scale_fill_gradient2(low = 'light blue', mid="yellow", high="red", midpoint=0.5) +
      labs(subtitle=paste0("Likelihood of introduction of ",focal_species," in 50 Years by drainage basin"), fill='Change in establishment risk') +
      xlab("Longitude") +
      ylab("Latitude")
  }

#---------------------------------------------------------------------#
# 7D. Internal validation plots
#
# We now get a couple of plots which show us uncertainty levels for
# the lakes which had introductions
#
#---------------------------------------------------------------------#

  if (plot_internal_validation == TRUE) {
    # Let's first plot the initial appearances
    Internal_Uncertainty <- ggplot(data=all_data_likelihoods[all_data_likelihoods$introduced == 1,]) + 
      geom_polygon(data=Norway1_sub_df, aes(long,lat,group=group), fill="grey") +
      geom_point(aes(x = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLongitude"], 
                     y = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLatitude"],
                     color=UncertaintyLevel)) + 
      scale_color_manual(values=c("black", "brown", "orange", "yellow")) +
      labs(subtitle=paste0("Current distribution of ",focal_species))
    
    
    
    Internal_Prediction <- ggplot(data=all_data_likelihoods[all_data_likelihoods$introduced == 1,]) + 
      geom_polygon(data=Norway1_sub_df, aes(long,lat,group=group), fill="grey") +
      geom_point(aes(x = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLongitude"], 
                     y = all_data_likelihoods[all_data_likelihoods$introduced == 1,"decimalLatitude"],
                     color=InitProbs_Class)) + 
      scale_color_manual(values=c("white", "yellow", "orange", "red")) +
      labs(subtitle=paste0("Current distribution of ",focal_species))
  }

maps <- list(predicted_appearances_pseudoAbs = predicted_appearances_pseudoAbs,
             predicted_appearances_fullAbs = predicted_appearances_fullAbs)

if (plot_internal_validation == TRUE) {
  maps[["Internal_Prediction"]] <- Internal_Prediction 
  maps[["Internal_Uncertainty"]] <- Internal_Uncertainty
  }

if (give_basin_plots == TRUE) {
  maps[["basin_map_pseudoAbs"]] <-basin_map_pseudoAbs
  maps[["basin_map_fullAbs"]] <- basin_map_fullAbs
}



saveRDS(maps, file=paste0("./Data/Species_Data/",gsub(' ','_',focal_species),"/maps.RDS"))



