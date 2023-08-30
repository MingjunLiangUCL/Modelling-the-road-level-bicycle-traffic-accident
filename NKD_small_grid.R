rm(list=ls())

# install devtools
install.packages('devtools')

# install the correct version of landscapemetrics (with the fixes I made)
devtools::install_github("r-spatialecology/landscapemetrics", force = TRUE, local = F)

library(sf)
library(spNetwork)
library(tmap)


setwd("/Users/liangmingjun/Desktop/UCL/Individual Project")
highways <- read.csv("merged_highways.csv")
highways <- highways[c("gml_id", "geometry", "length", "slope", 
                       "start_point_count", "end_point_count")]
highways$geometry <- st_as_sfc(highways$geometry)
highways_sf <- st_as_sf(highways)
#cycle_accident_2019_sf <- st_read("cycle_accident_2019_gdf.shp")
cycle_accident_2020_sf <- st_read("cycle_accident_2020_gdf.shp")



crs <- st_crs(27700)
highways_sf <- st_set_crs(highways_sf, crs)
#cycle_accident_2019_sf <- st_set_crs(cycle_accident_2019_sf, crs)
cycle_accident_2020_sf <- st_set_crs(cycle_accident_2020_sf, crs)


# Define the grid size
ny <- 9
nx <- 9
# Calculate the size of each cell
dy <- (205000-155000) / ny
dx <- (565000-500000) / nx
# Create a dataframe to store the results
results <- data.frame()

# Loop over the grid
for(i in 0:(ny-1)) {
  for(j in 0:(nx-1)) {
    if(i == 8 && (j == 4 || j == 5)) {
      next
    }
    
    cat(i, j, "\n")
    # Define the bounding box for this cell
    xmin <- 500000 + j*dx
    ymin <- 155000 + i*dy
    xmax <- xmin + dx
    ymax <- ymin + dy
    bbox <- st_sf(geometry = st_sfc(st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax), c(xmin, ymin)))), crs = 27700))
    
    # Clip the data to this bounding box
    clipped_highways <- st_intersection(highways_sf, bbox)
    if(nrow(clipped_highways) == 0) {
      next
    }
    # print(st_geometry_type(clipped_highways))
    clipped_highways <- clipped_highways[st_geometry_type(clipped_highways) == "LINESTRING", ]

    # clipped_cycle_accident_2019 <- st_intersection(cycle_accident_2019_sf, bbox)    
    clipped_cycle_accident_2020 <- st_intersection(cycle_accident_2020_sf, bbox)

    if (nrow(clipped_cycle_accident_2020) == 0) {
    # if (nrow(clipped_cycle_accident_2019) == 0) {
      next
    }

    # Perform the network kernel density calculation
    lixels <- lixelize_lines(clipped_highways,200,mindist = 50)
    samples <- lines_center(lixels)
  
    densities <- nkde(clipped_highways, 
                      events = clipped_cycle_accident_2020,
                      w = rep(1,nrow(clipped_cycle_accident_2020)),
                      samples = samples,
                      kernel_name = "quartic",
                      bw = 300, div= "bw", 
                      method = "discontinuous", digits = 1, tol = 1,
                      grid_shape = c(1,1), max_depth = 8,
                      agg = 5, #we aggregate events within a 5m radius (faster calculation)
                      sparse = TRUE,
                      verbose = FALSE)

    
    # Store the result
    if (any(is.na(densities))) {
      print("Found NA or NAN values in 'densities'")
    }
    samples$density <- densities
    results <- rbind(results, samples)
  }
}


print(results)


# rescaling to help the mapping
results$density <- results$density*1000
library(RColorBrewer)

results2 <- results[order(results$density),]
# write.csv(results2, file = "highways_nkd_centroid.csv", row.names = FALSE)
st_write(results2, "Data/highways_nkd_centroid/highways_nkd_centroid_2020.shp", append = FALSE)

print(length(results2))



colorRamp <- brewer.pal(n = 7, name = "Spectral")
colorRamp <- rev(colorRamp)

title <- paste0("bike accident density")
tm_shape(highways_sf) + 
  tm_lines("black") + 
  tm_shape(results2) + 
  tm_dots("density", style = "kmeans", palette = colorRamp, n = 7, size = 0.1) + 
  tm_layout(legend.outside = TRUE, 
            main.title = title , main.title.size = 1)


results_subset <- results2[, c("gml_id", "density")]
results_subset_df <- st_drop_geometry(results_subset)
merged_data <- merge(highways_sf, results_subset_df, by = "gml_id", all.x = TRUE)
merged_data$density[is.na(merged_data$density)] <- 0
merged_data <- st_zm(merged_data)
# write.csv(merged_data, file = "highways_nkd_line.csv", row.names = FALSE)
st_write(merged_data, "Data/highways_nkd_line/highways_nkd_line_2020.shp", append = FALSE)




# Zoom to a specific area
bbox <- st_as_sfc(st_bbox(c(xmin = 525000, ymin = 180000, xmax = 530000, ymax = 185000), crs = st_crs(highways_sf)))

highways_zoomed <- st_intersection(highways_sf, bbox)
results2_zoomed <- st_intersection(results2, bbox)

title <- paste0("bike accident density")
tm_shape(highways_zoomed) + 
  tm_lines("black") + 
  tm_shape(results2_zoomed) + 
  tm_dots("density", style = "kmeans", palette = colorRamp, n = 7, size = 0.1) + 
  tm_layout(legend.outside = TRUE, 
            main.title = title, 
            main.title.size = 1)


