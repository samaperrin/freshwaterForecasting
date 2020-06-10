#---------------------------------------------------------------------#
# 8. Consolidated visualisations ####
#
# Whilst script 7 allows us to look at distributions on a 
# species-by-species basis, this allows us to compare parameter effects
# across species.
# 
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 8A. Compare parameter effects
# 
# Show parameter effects by species, and by parameter
#---------------------------------------------------------------------#

paramInts_bySpecies <- list()
paramInts_byParam <- list()
allParams <- matrix(NA, nrow=0, ncol=5)

# First come up with lists and data frames showing all parameter effects by species
for (i in 1:length(species_list)) {
  all_analytics <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),
                                  "/whole_model_analytics_",population_threshold,"km_interaction.RDS"))
  each_species <- all_analytics$parameter_estimates
  each_species$species <- species_list[i]
  paramInts_bySpecies[[i]] <- each_species
  allParams <- rbind(allParams,each_species)
}
names(paramInts_bySpecies) <- species_list

allParams$commonNames <- rep(c("Whitefish", "Pike", "Roach", "Rudd", "Perch", "Rainbow trout"),
                            each= length(unique(allParams$variable)))

# for (j in 1:nrow(paramInts_bySpecies[[1]])) {
#   paramMat <- as.data.frame(matrix(NA, nrow = 6, ncol= 4))
#   colnames(paramMat) <- c("species", "lower", "mean", "upper")
#   for (i in 1:length(paramInts_bySpecies)) {
#     paramMat[i,] <- c(species_list[i],paramInts_bySpecies[[i]][j,3:5])
#   }
#   paramInts_byParam[[j]] <- paramMat
# }
# names(paramInts_byParam) <- paramInts_bySpecies[[1]][,2]



allParams$significance <- as.factor(ifelse(allParams$lower > 0 | allParams$upper < 0, 0,1))

# Need to give graphs better labels
# First create as many interactions as you need
for (i in 1:length(interaction_terms)) {
  int_names_one <- paste0("interaction",i)
  if(i == 1) {int_names <- int_names_one} 
  else {
    int_names <- c(int_names, int_names_one)
  }
}

if (!is.na(interaction_terms[[1]][1])) {
  supp.labs <- c(int_names,"A. Area", "B. Distance to Road","H. Downstream Populations","E. Human Footprint Index","F. Nearby Populations",
               "D. Population Distance","C. Temperature", "G. Upstream Populations")
  names(supp.labs) <- paramInts_bySpecies[[1]][,"variable"]
} else {
  supp.labs <- c("A. Area", "B. Distance to Road","H. Downstream Populations","E. Human Footprint Index","F. Nearby Populations",
                 "D. Population Distance","C. Temperature", "G. Upstream Populations")
  names(supp.labs) <- paramInts_bySpecies[[1]][,"variable"]
}

library(RColorBrewer)

# First figure, showing comparison of parameter effects
effects_byParam <- ggplot(allParams, aes(x=commonNames, y=mean,col=commonNames, linetype = significance)) +
  scale_colour_grey(aesthetics = "colour", start = 0, end = 0.8) + 
  ylim(-2,3) +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=0.25) +
  geom_point(aes(x=commonNames, y=mean, shape = commonNames)) +
  facet_wrap(~ variable, labeller = labeller(variable = supp.labs)) + 
  guides(linetype = FALSE) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right") + 
  geom_hline(yintercept=0)


# same figure, showing comparison of parameter effects with different grouping

effects_bySpecies <- ggplot(allParams, aes(x=variable, y=mean,col=variable, linetype = significance)) + 
  scale_color_brewer(palette = "Dark2") + 
  ylim(-1.2,1) +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=0.25) +
  facet_wrap(~ species) +
  guides(linetype = FALSE) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0)


#---------------------------------------------------------------------#
# 8B. Compare deviance values
# 
# Gets deviance values for our different models with our different
# values for radius in which to include nearby populations 
# (5, 10 or 20km).
#---------------------------------------------------------------------#


totalsMat <- as.data.frame(matrix(NA,nrow = length(species_list),ncol=4))
colnames(totalsMat) <- c("species","total","introduced","absent")

for (i in 1:length(species_list)) {
  all_data <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),"/whole_model_data_",population_threshold,"km.RDS"))
  totals <- all_data$raw_data %>% filter(native == 0) %>%
    group_by(introduced) %>% tally()
  totalsMat[i,] <- c(species_list[i],sum(totals$n),totals$n)
}  
  
devianceMat <- as.data.frame(matrix(NA,nrow = length(species_list),ncol=4))
colnames(devianceMat) <- c("species","5km","10km","20km")
for (i in 1:length(species_list)) {
  analysis_5km <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),"/whole_model_analytics_5km.RDS"))
  deviance5km <- analysis_5km$deviance
  analysis_10km <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),"/whole_model_analytics_10km.RDS"))
  deviance10km <- analysis_10km$deviance
  analysis_20km <- readRDS(paste0("./Data/Species_Data/",gsub(' ','_',species_list[i]),"/whole_model_analytics_20km.RDS"))
  deviance20km <- analysis_20km$deviance
  devianceMat[i,] <- c(species_list[i],deviance5km,deviance10km,deviance20km)
}  
devianceMat[, 2:4] <- sapply(devianceMat[, 2:4], as.numeric)
devianceMat[7,] <- c("Total", colSums(devianceMat[,2:4]))

#---------------------------------------------------------------------#
# 8C. Visualisation of lake density
# 
# Shows density of lakes throughout Norway
#
#---------------------------------------------------------------------#

total_lake_data <- readRDS(file=paste0("./Data/all_data.RDS"))

# Produce a map of Norway to use
Norway<-getData("GADM", country="NO", level=0)
Norway1<-getData("GADM", country="NO", level=1)
Norway1_sub<-Norway1[!(Norway1@data$NAME_1 %in% c('Troms', 'Finnmark', 'Nordland')),]
par(mar=c(1,1,1,1))
Norway_df <- fortify(Norway)
Norway1_sub_df <- fortify(Norway1_sub)
Norway1_sub_sf <- st_as_sf(Norway1_sub)

Norway_asOne <- st_union(Norway1_sub_sf)
#plot(Norway_asOne)

lake_density_Norway <- ggplot(data=total_lake_data) +
  geom_sf(data = Norway_asOne, fill = "white") +
  geom_hex(aes(x = total_lake_data[,"decimalLongitude"], 
               y = total_lake_data[,"decimalLatitude"]),
           bins = n_bins) +
  scale_fill_gradient(low="white", high="dark blue", #colors in the scale
                       #midpoint=40,    #same midpoint for plots (mean of the range)
                       breaks=seq(0,40,40/4), #breaks in the scale bar
                       limits=c(0,40)) + 
  labs(fill="Lakes") +
  xlab("Longitude") +
  ylab("Latitude")
# 
# Sweden<-getData("GADM", country="SE", level=0)
# Sweden1<-getData("GADM", country="SE", level=1)
# Sweden1_sub<-Sweden1[!(Sweden1@data$NAME_1 %in% c('Västerbotten', 'Norrbotten')),]
# par(mar=c(1,1,1,1))
# Sweden_df <- fortify(Sweden)
# Sweden1_sub_df <- fortify(Sweden1_sub)
# Sweden1_sub_sf <- st_as_sf(Sweden1_sub)
# Sweden1_sub_sf <- st_transform(Sweden1_sub_sf,25833)
# 
# basin_geom_sf <- readRDS("./Data/waterRegions.rds")
# 
# # And plot
# basin_map_noIncrease <- ggplot() + 
#   geom_sf(data = Norway1_sub_sf,fill="light green", lwd=0) +
#   geom_sf(data = Sweden1_sub_sf,fill="light blue", lwd=0) +
#   geom_sf(data = basin_geom_sf,fill=NA) +


comparative_stats <- list(parameter_effects = allParams, effectsFigure1 = effects_bySpecies, effectsFigure2 = effects_byParam,
                          deviance_comparison = devianceMat, lakeDensity = lake_density_Norway)
saveRDS(comparative_stats,file="./Data/ComparativeStats.RDS")
  
  
  
