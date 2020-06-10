### SUPPLEMENTARY INFORMATION FOR PERRIN ET AL., 2020


### This is the overarching script that ties everything together. You can of course
### dive into the source scripts, but theoretically you can run all the data gathering 
### and analysis from this script.

### Let's get all the libraries in
packages.needed <- c('dplyr', 'pool', 'postGIStools', 'sf', 'getPass', 
                     'rgbif', 'mapedit', 'rio', 'raster', 'FNN', 'readr',
                     'ggplot2', 'greta', 'stringr', 'lwgeom', 'purrr')

packages.toinstall <- setdiff(packages.needed, rownames(installed.packages())  
) 
if(length(packages.toinstall)) install.packages(packages.toinstall) 
lapply(packages.needed, library, character.only = TRUE)

### I also have a few custom-made functions that you'll need
source("ExtraFunctions.R")


### Script 1 downloads occurrence data for all your species and matches it to the nearest lake. 

# Define species of interest

species_list <- c("Coregonus lavaretus", "Esox lucius", "Rutilus rutilus",
                  "Scardinius erythrophthalmus", "Perca fluviatilis", "Oncorhynchus mykiss")
common_names <- c("Whitefish", "Pike", "Roach",
                  "Rudd", "Perch", "Rainbow trout")

# All the data you create will be transferred into a folder. Let's create that folder now.
if (dir.exists(paste0("./Data")) == FALSE
) {dir.create(paste0("./Data"))}

# The first download from GBIF is tricky. You need to intiate a download, which can be 
# done using script 1. As the script should spit out though, this can take a while.
# So you'll have to wait a bit after running the script for the first time. After you're done waiting,
# set 'initiate_download' to FALSE.

initiate_download <- FALSE            # Have you already initiated a download? If so, set to false.
GBIF_download <- FALSE                # Have you already downloaded occurrence data from GBIF? If so, set to false.
download_lakes <- FALSE               # Have you already downloaded the lakes? If so, set to false.
get_all_lakes <-  FALSE               # Have you downloaded all lake data from NOFA? If so, set to false.
delete_north <-  FALSE                # Do you want to use lakes north of Tronderlag? If so, set to true.
download_native_range <- FALSE        # Have you already downloaded species' native ranges? If so, set to false.

# Data from GBIF has a geographic coordinate which is then matched to the nearest lake.
# dist_threshold sets the distance between a point and the nearestlake which
# is an acceptable match. So if dist_threshold is set to 50, a point will only be
# matched to a lake if it is less than 50m away.
dist_threshold <- 50

source("./R/1_GetIntroductions.R")

# Next step is a fairly short script to add in all environmental covariates which can be 
# calculated regardless of species. Things like lake area, surface temperature,
# Human Footprint Index.

size_threshold <- 0.02               # All lakes below this size in km2 will be discarded.

HFP_download <- FALSE                # Set this to false if you already have raster data
                                     # showing the Human Footprint Index downloaded in your data folder.
catchment_download <- FALSE          # Set this to false if you already have the lakes sorted by catchment and the 
                                     # catchment geometries downloaded. This takes AGES, so the default is always false here.

source("./R/2_BioticDataAddition.R")

# From now on everything becomes species specific. The next script gets connectivity data 
# for indivudal species.
focal_species <- species_list[6]

# Define radius you want to measure nearby populations by (see main text).
population_threshold <- 20

# Don't worry about the duplication if the only lake that has been duplicated is 39447.
# If there are others, let me know.

source("./R/3_SpeciesDataAddition.R")

# Now we run the preliminary model. Only thing we need to choose is what our size limit 
# on lakes will be.

# Parameters used in the model are as follows.
# 1. Lake_area
# 2. Distance to nearest road
# 3. Average temperature of warmest quarter
# 4. Euclidean distance to the nearest population
# 5. Human footprint index
# 6. Number of populations within km radius specified by 'population threshold'
# 7. Number of populations upstream
# 8. Number of populations downstream

# To leave one or more of these parameters out, simply input the number into the object below.
parameters_to_ignore <- NA
interaction_terms <- list(c(5,6), c(5,4))

# This script will take a while, as it's running a model on up to 40,000 lakes. Grab a 
# coffee. Teach it to do algebra. There are two options, one which includes an interaction effect
# and one which doesn't. That which does is our default.
#source("./R/4_FullModelConstruct.R")
source("./R/4_FullModelConstruct_InteractionsBeta.R")

# Next one gives you uncertainty from each lake, convergence diagnostic,
# beta intervals, and model deviance. You can check them out individually by
# looking at the object 'model_analysis' after they've run.

#source("./R/5_ModelAnalysis.R")
source("./R/5_ModelAnalysis_InteractionsBeta.R")

# Now we can run our simulations. need to define which model we want to use.
# n_loops simply tells us how many iterations we want to run. Be careful,
# because the run time can be enormous. With 40,000 lakes, 100 loops generally
# takes 10 minutes. The temp_scenario object dictates whether or not 
# you want to introduced an increase in temperature to the scenario.
# We advise against doing this.
n_loops <- 1000
temp_scenario <- FALSE

#source("./R/6_Looping.R")
source("./R/6_Looping_InteractionsBeta.R")

# The following scripts are slight variations on the script above, in that
# they do our retrospective scenario modelling.

#source("./R/6a_HistoricalForecast.R")
source("./R/6a_HistoricalForecast_InteractionsBeta.R")

# Last 2 scripts are all about data visualisation. n_bins dictates resolution of
# your hexagons. If you have already downloaded basin data, let download_basins =
# FALSE.

n_bins <- 120
download_basins=FALSE

# The script produces an object calls 'maps', which you can then go through to find
# the maps you're after. These also include maps which are based on watersheds,
# as well as a few which reflect internal validation.

source("./R/7_Visualisations.R")

# Lastly, we have some comparative visualisations, including beta parameter effects
# and deviance comparisons. To take a look after the script has run, check the object
# comparative_stats.
# source("./R/8_ComparativeViz.R")
source("./R/8_ComparativeViz_InteractionsBeta.R")


