extract_raster <- function(path = getwd(), pattern = ".tif$", coords, check.packages = T, select = "all", wsize = 5) {
  if (check.packages) { 
    if (!require("pacman")) install.packages("pacman") 
    pkgs = c("raster", "sp", "pbapply")
    pacman::p_load(char = pkgs, character.only = T)
  } 
  raster_list <- list.files(path = path, recursive = T, pattern = pattern, full.names = TRUE)
  col_names <- tools::file_path_sans_ext(basename(raster_list))
  if(select[1] != "all") raster_list <- raster_list[which(col_names %in% select)]
  col_names <- tools::file_path_sans_ext(basename(raster_list))
  
  # Function to fill NA
  #fill.na <- function(x, i =  ceiling(wsize^2 / 2)) {
  #  if( is.na(x)[i] ) {
  #    return( round(mean(x, na.rm=TRUE),0) )
  #  } else {
  #    return( round(x[i],0) )
  #  }
  #}  
  
  df <- pbsapply(raster_list, FUN = function(x) {
    temp <- raster(x) 
    if (proj4string(temp) != proj4string(coords)) coords <- spTransform(coords, CRSobj = crs(temp))
    #temp <- crop(temp, coords)
    #temp <- focal(temp, w = matrix(1, wsize, wsize), fun = fill.na, pad = TRUE, na.rm = FALSE)
    df <- raster::extract(temp, coords)
    })
  df <- as.data.frame(df)
  names(df) <- col_names
  return(df)
} 

