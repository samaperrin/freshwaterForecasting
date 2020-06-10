#---------------------------------------------------------------------#
# 1. Match GBIF Data to Lakes and Native Ranges
#---------------------------------------------------------------------#

# This script downloads occurrence data from GBIF, matches it to a suitable lake, then decides whether
# or not the occurrences are introductions based on whether or not they occur in the species' native range.


#---------------------------------------------------------------------#
# 1A. Get data from GBIF ####
#
# Load occurrences from GBIF. Initiates download in GBIF, then
# pulls everything into the R environment.
# 
#---------------------------------------------------------------------#

# THe following initiates a download of your data. If it hasn't already been
# downloaded, it will require you to run the script again later when the
# download has finished.

if (initiate_download == TRUE) {
  # set GBIF user account credentials 
  options(gbif_user=rstudioapi::askForPassword("my gbif username"))
  options(gbif_email=rstudioapi::askForPassword("my registred gbif e-mail"))
  options(gbif_pwd=rstudioapi::askForPassword("my gbif password"))
  
  keys <- as.integer()
  for(i in 1:length(species_list)){
    keys[i] <- name_suggest(q=species_list[i], rank='species')$key[1] 
  }
 
  
  # Take the key given by the code above and initiate the download using the following command
  download_key <- occ_download(
    paste0('taxonKey = ',paste(keys[1:6],collapse=",")),
    type = "or"
  ) %>% 
    occ_download_meta
  
  # store download_key for re-use 
  saveRDS(download_key,"./Data/download_key.rds")
  
  stop("Download being initiated. You will have to wait for download to complete in GBIF. Suggest that
       you wait 15 minutes then run script again with initiate_download set to FALSE.")
  
} else {download_key <- readRDS("./Data/download_key.rds")}

# Now, download lakes into Data folder

if (GBIF_download == TRUE) {
  temp <- tempdir()
  download.file(url=download_key$downloadLink,
                destfile=paste0(temp,"/tmp.zip"),
                quiet=TRUE, mode="wb")
  
  # Unzip dowloaded occurrence file. We also need to narrow down to Norwegian data and get rid of
  # all occurrences marked with 'absence', and take only relevant columns.
  species_distribution <- rio::import(unzip(paste0(temp,"/tmp.zip"),files="occurrence.txt"))
  species_distribution_all <- species_distribution %>%
    filter(countryCode == "NO" & occurrenceStatus != "absent") %>%
    dplyr::select(gbifID,scientificName, occurrenceID, occurrenceStatus, catalogNumber,decimalLongitude, decimalLatitude,species,taxonKey,datasetKey, locality,municipality,county,countryCode,locationID,
                  eventDate,year,month,day,samplingProtocol,eventID,fieldNumber,
                  recordedBy,dynamicProperties,collectionCode,datasetName,license,institutionCode)
  
  saveRDS(species_distribution_all,"./Data/allSpecies_GBIFDownload.RDS")
  unlink(temp)
} else {species_distribution_all <- readRDS("./Data/allSpecies_GBIFDownload.RDS")}


print("Finished importing raw species data")

#---------------------------------------------------------------------#
# 1B. Download lake data and match to occurrences ####
#
# Load a series of lake polygons and match occurrences to nearest 
# lakes.
#---------------------------------------------------------------------#


## Get occurrence data in spatial data format

# Need to match the point data to the occurrence data now
species_dist_short <- species_distribution_all[complete.cases(species_distribution_all$decimalLongitude),]

species_dist_short$longitude <- species_dist_short$decimalLongitude
species_dist_short$latitude <- species_dist_short$decimalLatitude
species_dist_short$scientificNameShort <- word(species_dist_short$scientificName, 1,2, sep=" ")

# Have to convert it using the sf package first. Standard CRS uses the EPSG 32633.
dist_sf <- st_as_sf(species_dist_short, coords = c("longitude", "latitude"), 
                    crs = 4326)
dist_sf <- st_transform(dist_sf, 32633)


# Following code brings the lake map down off box, so you'll need access to the box folder.
if (download_lakes == TRUE) {
  temp <- tempdir()
  download.file("https://ntnu.box.com/shared/static/6vu4de2birf9onmej2gorexwaxposuaa.zip",destfile = paste0(temp,"/NVE_innsjodatabasen.zip"))
  unzip(paste0(temp,"/NVE_innsjodatabasen.zip"),exdir=temp)
  lakes <- st_read(paste0(temp,"/Innsjo_Innsjo.shp"))
  saveRDS(lakes,"./Data/lake_polygons.rds")
} else {
  lakes <- readRDS("./Data/lake_polygons.rds")
}


# Now we join 2 columns to our data frame which give the number of the nearest lake and 
# the distance to that lake.
occ_with_lakes <- st_join(dist_sf, lakes, join = st_nearest_feature)
index <- st_nearest_feature(x = occ_with_lakes, y = lakes) # index of closest lake
closest_lakes <- lakes %>% dplyr::slice(index) # slice based on the index
dist_to_lake <- st_distance(x = occ_with_lakes, y= closest_lakes, by_element = TRUE) # get distance
occ_with_lakes$dist_to_lake <- as.numeric(dist_to_lake) # add the distance calculations to match data

# Get rid of all coords which are more than a certain distance from the nearest lake
occ_OKlakes <- occ_with_lakes %>% filter(dist_to_lake < dist_threshold)

saveRDS(occ_OKlakes,"Data/allSpecies_GBIFDownload_Edited.RDS")

print("Finished matching species distribution to native range and lakes")


#---------------------------------------------------------------------#
# 1C. Pull other lakes from NOFA to establish absences ####
#
# Loads all other lakes from NOFA so that we have absences to match 
# against presences eventually.
#---------------------------------------------------------------------#

if (get_all_lakes ==  TRUE) {
  
  pg_user=rstudioapi::askForPassword("Wallace username")
  pg_password=rstudioapi::askForPassword("password")
  
  # Connect to nofa, extract polygons
  pool <- dbPool(drv = RPostgreSQL::PostgreSQL(), dbname = 'nofa',
                 host = 'vm-srv-wallace.vm.ntnu.no', user = pg_user, password = pg_password,
                 idleTimeout = 36000000,
                 options="-c search_path=nofa"
  )
  con <- poolCheckout(pool)
  
  # Get lakes for everything
  all_lakes <- tbl(con,"location") %>%
    filter(countryCode =='NO') %>%
    dplyr::select(locationID, waterBodyID, county, no_vatn_lnr, decimalLatitude, decimalLongitude) %>%
    collect()
  saveRDS(all_lakes,file="./Data/all_lakes.RDS")
} else {
  all_lakes <- readRDS("./Data/all_lakes.RDS")
}

print("Downloaded all lakes.")

#---------------------------------------------------------------------#
# 1D. Download native distribution ####
#
# Load native distribution range and assign establishmentMeans status
# to occurrences accordingly 
#---------------------------------------------------------------------#

## Get native distribution from https://doi.org/10.21400/1mwt3950 (NB: data are in EPSG:32633)
if (download_native_range == TRUE) {
  
  url <- "https://api.loke.aws.unit.no/dlr-gui-backend-resources-content/v2/contents/links/1efa6a5a-74b7-46d1-9bd6-49a7b3c58d030e379585-05d5-42ed-b00d-35d2d2b0f9bf1fc5a283-f739-4205-8c84-9748cb6d76ae"
  temp <- tempdir()
  download.file(url,paste0(temp,"/hk_native.zip"))
  unzip(paste0(temp,"/hk_native.zip"),exdir=temp, overwrite=TRUE)
  files <- unzip(paste0(temp,"/hk_native.zip"),list = TRUE)[,1]
  shapefile <- files[str_detect(files,".shp")]
  hk_distribution_map <-  st_read(paste0(temp,"/",shapefile))
  hk_distribution_map <- st_transform(hk_distribution_map, 4326) # reproject data to lat/long wgs84
  
  saveRDS(hk_distribution_map,"Data/native_distribution.rds")
  
} else {
  hk_distribution_map <- readRDS("Data/native_distribution.rds")
}


# Figure out which lakes are in the native range or not (this might take a while)

# Dulicate coordinates for later so that we can have a geometry column (produced when
# using the sf package) and still have lat and long columns.
all_lakes$latitude <- all_lakes$decimalLatitude
all_lakes$longitude <- all_lakes$decimalLongitude

# Convert to sf file
all_lakes_sf <- st_as_sf(all_lakes, coords = c("longitude", "latitude"), 
                         crs = 4326)
all_lakes_sf <- st_transform(all_lakes_sf, 32633)

# See which lakes are in which native range. Need to create a special case for Rainbow Trout, as
# it doesn't have a native range in Norway

for(i in 1:length(species_list)) {
  if (species_list[i] != 'Oncorhynchus mykiss') {
    fish_dist <- hk_distribution_map[(hk_distribution_map$canonicaln==species_list[i]
                                      & hk_distribution_map$establishm=="native"),]
    fish_dist <- st_transform(fish_dist, 32633)
    inOut <- st_intersects(all_lakes_sf, fish_dist)
    
    # This gives us a list, which gives the status of occurrence for each lake. Since we want
    # no occurrence status, the list elements with length 0 are the ones that are not in the 
    # native range.
    all_lakes_sf[,paste0(gsub(" ","_",species_list[i]),"_native")] <- lengths(inOut)
  } else {
    all_lakes_sf[,paste0(gsub(" ","_",species_list[i]),"_native")] <- 0
  }
}

# Need to create a bunch of columns showing presence and introductions for all species. Function runs 
# much quicker when converted to data frames instead of sf objects.
indicator <- as.data.frame(occ_OKlakes)
all_lakes_df <- as.data.frame(all_lakes_sf)

for (i in 1:length(species_list)) {
  present_lakes <- indicator[indicator$scientificNameShort==species_list[i],"vatnLnr"]
  all_lakes_sf[,paste0(gsub(" ","_",species_list[i]),"_presence")] <- ifelse(all_lakes_df$no_vatn_lnr %in% present_lakes, 1, 0)
  all_lakes_sf[,paste0(gsub(" ","_",species_list[i]),"_introduced")] <- 
    ifelse(all_lakes_df$no_vatn_lnr %in% present_lakes & 
             all_lakes_df[,paste0(gsub(" ","_",species_list[i]),"_native")] == 0, 1, 0)
}



# Lastly, Get rid of all lakes north of Trondelag
if (delete_north == TRUE) {
  all_lakes_sf <- all_lakes_sf %>%
    filter(!(county %in% c("Troms","Finnmark","Nordland")))
}
print("Finished editing data")

saveRDS(all_lakes_sf, file=paste0("./Data/introductions.RDS"))

