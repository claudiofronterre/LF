# Load required packages and functions
if (!require("pacman")) install.packages("pacman")
pkgs = c("sp", "mapview", "tmap") # package names
pacman::p_load(pkgs, character.only = T)

# Load and clean data --------------
lf <- rgdal::readOGR("data/LF_SiteLevel_Dataset_Feb2018.shp")

# Change projection (use one tha preserves the distance)
proj4string(lf)
lf <- spTransform(lf, CRS("+init=epsg:3857 +units=km"))

# Visualisation -----------------

# Interactive visualisation of data
mapView(lf, zcol = c("PoT", "Prevalence", "Method_1"), legend = F, burst = T, hide = T)

# Static visualization of data
africa <- rgdal::readOGR("data/Africa.shp")
proj4string(africa) <- CRS("+init=epsg:4326") 
africa <- spTransform(africa, proj4string(lf))

map <- tm_shape(africa) +
  tm_borders("black") +
  tm_fill("white") +
       tm_shape(lf) +
  tm_symbols(col = "Method_1", title.col = "Survey data by\ndiagnostic method", 
             size = .05, border.col = "black", palette = c("red", "blue")) +
  tm_compass(position = "left") +
  tm_layout(bg.color = "lightblue", design.mode = F, legend.bg.color = "white", legend.frame = "black") +
  tm_scale_bar() 
map      

