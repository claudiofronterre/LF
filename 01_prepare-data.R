# Load required packages and functions
if (!require("pacman")) install.packages("pacman")
pkgs = c("sp", "mapview", "tmap", "raster") # package names
pacman::p_load(pkgs, character.only = T)

# Load and clean data --------------
lf <- read.csv("data/LF_Site_Level_Data_Feb2018.csv")

# Remove inconsistent and not useful points
lf <- lf[!is.na(lf$Examined), ] # Remove sites with Examined = NA
lf <- lf[lf$Examined >= lf$Positive, ] # Remove sites with Examined < Positive
lf <- lf[lf$Examined != 0, ] # Remove sites with Examined = 0
lf <- lf[lf$Method_0 == "Mapping survey" & lf$PoT == "Preintervention", ]

# Check if important columns have NA and remove them
View(apply(lf, 2, function(x) sum(is.na(x))))

#

# Change projection (use one that preserves the distance)
lf_sp <- lf
coordinates(lf_sp) <- ~ Longitude + Latitude
proj4string(lf_sp) <- CRS("+init=epsg:4326")
lf_sp <- spTransform(lf_sp, CRS("+init=epsg:3857 +units=km"))

# Save the clean file 
#saveRDS(as.data.frame(lf_sp), file = "data/lf_clean.rds")

# Visualisation -----------------

# Interactive visualisation of data
mapView(lf_sp, zcol = "Examined", legend = T, at = c(1, 50, 100, 500, 1000, 5000), cex = "Examined")
mapView(lf_sp, zcol = c("PoT", "Prevalence", "Method_1"), legend = T, burst = T, hide = T)

# Static visualization of data
africa <- rgdal::readOGR("data/Africa.shp")
proj4string(africa) <- CRS("+init=epsg:4326") 
africa <- spTransform(africa, proj4string(lf_sp))

map <- tm_shape(africa) +
  tm_borders("black") +
  tm_fill("white") +
       tm_shape(lf_sp) +
  tm_symbols(col = "Method_1", title.col = "Survey data by\ndiagnostic method", 
             size = .25, border.col = "black", palette = c("red", "blue")) +
  tm_compass(position = "left") +
  tm_layout(bg.color = "lightblue", design.mode = F, legend.bg.color = "white", legend.frame = "black") +
  tm_scale_bar() 
map      

# tmap_mode(mode = "view")
# map
# ttm()

save_tmap(tm = map, filename = "figs/map.pdf", width = 10, height = 12)

# Select points where both diagnostics where used ------------------------------------------------
ict <- which(lf$Method_1 == "Serological")
mf <- which(lf$Method_1 == "Parasitological")
lf_ict <- lf[ict, ]
lf_mf <- lf[mf, ]
lf_both <- merge(x = lf_mf, y = lf_ict, by = c("Country", "Latitude", "Longitude", "SurveyYear", "Examined"), 
                 suffixes = c("_mf", "_ict"))
lf_both$Country <- factor(lf_both$Country)
table(lf_both$Country)
lf_both <- lf_both[, c("Country", "Latitude", "Longitude", "SurveyYear", "Examined", "Positive_ict", "Positive_mf")]
lf_both$Prevalence_ict <- lf_both$Positive_ict/lf_both$Examined
lf_both$Prevalence_mf <- lf_both$Positive_mf/lf_both$Examined
saveRDS(lf_both, file = "data/lf_both.rds")

# Prepare prediction  grid
map <- raster("data/Mfprev_filled.tif")
map10km <- aggregate(map, fact = 10, filename = "data/map10km.tif", format = "GTiff", overwrite=TRUE)
map5km <- disaggregate(map10km, fact = 2, filename = "data/map5km.tif", format = "GTiff", overwrite=TRUE)
pred_coords <- rasterToPoints(map5km)[,c(1, 2)]/1000
saveRDS(pred_coords, file = "data/pred_coords.rds")

# Prepare covariates from raster (path is the location of all the rasters) 
coords <- data.frame(x = lf$Longitude, y = lf$Latitude) # extraction points
coordinates(coords) <- ~ x + y
proj4string(coords) <- CRS("+init=epsg:4326")
coords <- spTransform(coords, CRS("+init=epsg:3857"))
source("R/extract_raster.R")
covariates <- extract_raster(path = "~/Documents/Raster_Covariates/", coords = coords)
saveRDS(covariates, file = "data/covariates.rds")
