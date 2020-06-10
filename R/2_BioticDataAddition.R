#---------------------------------------------------------------------#
# 2. Import Biotic Data Relevant to All Species
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# 2A. Import Biotic Data ####
#
# At the moment this downloads lake area, perimeter, distance to road
# and temeprature.
# 
#---------------------------------------------------------------------#

# This script downloads biotic data that is common to all lakes.

# Import data we have so far
all_lakes_sf <- readRDS(paste0("./Data/introductions.RDS"))

pg_user=rstudioapi::askForPassword("Wallace username")
pg_password=rstudioapi::askForPassword("password")

# Connect to the NOFA database
pool <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
               host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
               idleTimeout = 36000000,
               options="-c search_path=nofa"
)
con <- poolCheckout(pool)

# And one more to get catchment data
pool_catchment <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
                         host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
                         idleTimeout = 36000000,
                         options="-c search_path=public,catchments"
)
con_catchment <- poolCheckout(pool_catchment)

# And another to get temperature data
temp_pool <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
                    host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
                    idleTimeout = 36000000,
                    options="-c search_path=environmental"
)
temp_con <- poolCheckout(temp_pool)

# Get area, shoreline, distance to road
biotic_data_import <- tbl(con, "lake") %>%
  dplyr::select(id,area_km2, perimeter_m, distance_to_road) %>%
  filter(id %in% !!all_lakes_sf$waterBodyID) %>%
  collect()

# Get temperature
temp_data_import <- tbl(temp_con, "lake_EuroLST_BioClim") %>%
  dplyr::select(waterBodyID, eurolst_bio10) %>%
  filter(waterBodyID %in% !!all_lakes_sf$waterBodyID) %>%
  collect()

# And get location data as well
waterRegion_import <- tbl(con, "location") %>%
  dplyr::select(waterBodyID, eb_waterregionID, ebint, maximumElevationInMeters) %>%
  filter(waterBodyID %in% !!all_lakes_sf$waterBodyID) %>%
  collect()
waterRegion_import <- waterRegion_import[!duplicated(waterRegion_import),]

print("Downloaded biotic data, starting on downloading geometries.")

# Import catchment size for use in getting rid of duplicates
catchment_import <- tbl(con_catchment,"lake_catchments") %>%
  dplyr::select(ebint, geom, catchment_area_km2) %>%
  filter(ebint %in% !!waterRegion_import$ebint) %>%
  collect()

print("Downloaded catchment geometries.")

#---------------------------------------------------------------------#
# 2B. Import HFP Data ####
#
# I've uplaoded the HFP (Human Footprint Index) data to the box folder 
# in raster form. Might need a better solution in the future. For more
# info on the HFP index click the link below.
# https://wcshumanfootprint.org/
# 
#---------------------------------------------------------------------#

# Get HFP data from box download

if (dir.exists(paste0("./Data/HFP")) == FALSE
) {dir.create(paste0("./Data/HFP"))}

if (HFP_download == TRUE) {
  temp <- tempdir()
  download.file("https://ntnu.box.com/shared/static/prw0tsdp9hucsyrkm3lml4oospx4f01e.zip",destfile = paste0(temp,"/HFP1993.zip"))
  unzip(paste0(temp,"/HFP1993.zip"),exdir="./Data/HFP")
  HFP_raster <- raster("./Data/HFP/HFP1993.tif")
  saveRDS(HFP_raster,"./Data/HFP_raster.rds")
} else {
  HFP_raster <- readRDS("./Data/HFP_raster.rds")
}

# Need to get site points in a format we can compare to raster data
site_points <- as.data.frame(all_lakes_sf)[,c("decimalLatitude","decimalLongitude")]
coordinates(site_points) = ~ decimalLongitude+ decimalLatitude
crs(site_points) <-  "+init=epsg:4326 +proj=longlat"

site_points_conv <- spTransform(site_points, CRSobj = crs(HFP_raster))

HFPValues <- extract(HFP_raster, site_points_conv)

print("Downloaded HFP data")

#---------------------------------------------------------------------#
# 2B. Build consolidated data frame ####
#
# This brings all the above data together and gets rid of any duplicates
# we might have.
# 
#---------------------------------------------------------------------#


all_lakes_dataBuild <- all_lakes_sf
all_lakes_dataBuild$HFP <- HFPValues                     # Import HFP data

# Area, shoreline, distance to road
all_lakes_dataBuild <- merge(all_lakes_dataBuild, biotic_data_import,by.x = "waterBodyID", by.y= "id", all.x=TRUE)  

# At this point we can speed things up by getting rid of lakes under our threshold for size
all_lakes_dataBuild <- all_lakes_dataBuild %>%
  filter(area_km2 > size_threshold)

# Aaaand temperature
all_lakes_dataBuild <- merge(all_lakes_dataBuild, temp_data_import, by = "waterBodyID", all.x=TRUE)

# And location data
all_lakes_dataBuild <- merge(all_lakes_dataBuild, waterRegion_import, by = "waterBodyID", all.x=TRUE)

# And catchment area
all_lakes_dataBuild <- merge(all_lakes_dataBuild, catchment_import[,c("ebint","catchment_area_km2")], by = "ebint", all.x=TRUE)

# Create utm values for later use
all_lakes_utm <- unlist(all_lakes_dataBuild$geometry) %>%
  matrix(ncol=2,byrow=TRUE) %>% as.data.frame()

all_lakes_dataBuild$utm_x <- all_lakes_utm[,1]
all_lakes_dataBuild$utm_y <- all_lakes_utm[,2]

# Get rid of geometry column.
all_lakes_final <- all_lakes_dataBuild %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  arrange(waterBodyID)

# Get only complete cases
all_lakes_final_complete <- all_lakes_final[complete.cases(all_lakes_final),]

# Get rid of all duplicates
all_lakes_final_complete <- all_lakes_final_complete[!duplicated(all_lakes_final_complete),]

# Some water bodies have different catchments for some reason. We take the bigger catchment
df.agg <- aggregate(catchment_area_km2 ~ waterBodyID, all_lakes_final_complete, max)
all_lakes_final_complete <- merge(df.agg, all_lakes_final_complete, all.x=TRUE, by = c("catchment_area_km2","waterBodyID"))

# Need to sort out our duplicated lakes
duplicated_lakes_vatnLnr <- all_lakes_final_complete[duplicated(all_lakes_final_complete$no_vatn_lnr),"no_vatn_lnr"]
lake_dupes <- all_lakes_final_complete[all_lakes_final_complete$no_vatn_lnr %in% duplicated_lakes_vatnLnr,]


# We have 4 lakes which have been duplicated, for the following reasons:
# 39447 - 2 different lakes, can just keep them both in there
# 41167 - Same lake, just delete one
lakeRow_41167 <- which(all_lakes_final_complete$no_vatn_lnr==41167)[2]
all_lakes_final_complete <- all_lakes_final_complete[-lakeRow_41167,]
# 24083 - Same lake spread across county border, merge HFP and you're fine
HFP_lake24083 <- mean(all_lakes_final_complete[all_lakes_final_complete$no_vatn_lnr == 24083, "HFP"])
all_lakes_final_complete[all_lakes_final_complete$no_vatn_lnr == 24083, "HFP"] <- HFP_lake24083
lakeRow_24083 <- which(all_lakes_final_complete$no_vatn_lnr==24083)[-1]
all_lakes_final_complete <- all_lakes_final_complete[-lakeRow_24083,]
# 7 - This is the problem one, the one with area of 19.6 is the one that should have presences
lakeRow_7 <- which(all_lakes_final_complete$no_vatn_lnr==7 & all_lakes_final_complete$area_km2 < 19)
all_lakes_final_complete <- all_lakes_final_complete[-lakeRow_7,]

# Get rid of the rest for now
duplicated_lakes_vatnLnr_removal <- all_lakes_final_complete[duplicated(all_lakes_final_complete$no_vatn_lnr),"no_vatn_lnr"]
rows_toGo <- which(all_lakes_final_complete$no_vatn_lnr %in% duplicated_lakes_vatnLnr_removal)
all_lakes_final_complete <- all_lakes_final_complete[-rows_toGo,]


saveRDS(all_lakes_final_complete, file=paste0("./Data/all_data.RDS"))

print("Data saved and can be found in Data folder (filename is all_data). Running catchment geometries now.")

#---------------------------------------------------------------------#
# 2D. Build geoms for water basin and catchment ####
#
# Need to build geometries for all water basins and catchments for later use
# 
#---------------------------------------------------------------------#

if (catchment_download == TRUE) {

# Have to transform this BACK into a geometry
all_lakes_forSpatial <- st_as_sf(all_lakes_final_complete, coords = c("decimalLongitude", "decimalLatitude"), 
                                                 crs = 4326)
all_lakes_forSpatial <- st_transform(all_lakes_forSpatial, 32633)

# Create vector of catchment characters to be converted into geoms
catchment_list <- vector()
for(i in seq_len(nrow(catchment_import))) {
  wkb_int <- structure(list(catchment_import$geom[i]), class = "WKB")
  catchment_list <- c(catchment_list,wkb_int)
}
# Conversion
catchment_geom <- st_as_sfc(catchment_list, EWKB = TRUE)
catchment_geom_sf <- catchment_geom %>%
  st_sf %>%
  st_cast
catchment_geom_sf$ebint <- catchment_import$ebint

# Now let's narrow down to only include the basins we need
catchment_geom_sf <- catchment_geom_sf %>%
  filter(ebint %in% unique(all_lakes_forSpatial$ebint))

# Align crs with our lakes
catchment_geom_sf <- st_transform(catchment_geom_sf, 32633)

# Create list
catchment_list2 <- list()

# Alot lakes to catchments they are in
time1 <- Sys.time()
for (i in seq_len(nrow(catchment_geom_sf))) {
  inOut <- st_intersects(all_lakes_forSpatial, catchment_geom_sf[i,])
  rows_toUse <- which(lengths(inOut) != 0)
  catchment_list2[[i]] <- as.data.frame(all_lakes_forSpatial)[rows_toUse,"waterBodyID"]
  time2 <- Sys.time()
  difftime_1 <- round(as.numeric(difftime(time2, time1,
                                          units = "mins")),2)
  if (i %% 100 == 0) {print(paste0("Run ", i, " finished in ",difftime_1, " minutes. Estimated ", 
                                 round(difftime_1*(nrow(catchment_geom_sf))/i-difftime_1,2), " minutes left.") )}
  }

# Give catchments names
names(catchment_list2) <- catchment_geom_sf$ebint

catchment_data <- list(catchmentsByLake = catchment_list2, catchmentGeometries = catchment_geom_sf)


saveRDS(catchment_data,file="./Data/lakesInCatchments.RDS")
}
