rm(list=ls())

library(sf)
library(CARBayes)



setwd("/Users/liangmingjun/Desktop/UCL/Individual Project")
# highways_centroid <- read.csv("highways_nkd_centroid.csv", row.names=NULL)
# file_path <- "Data/highways_nkd_centroid/highways_nkd_centroid.shp"
file_path <- "Data/final_highways/final_highways_2020.shp"
highways_2020 <- st_read(file_path)
names(highways_2020)
highways_2020 <- highways_2020[c("gml_id", "slope", "length", "exit_1", "exit_2", 
                                 "exit_3", "exit_4", "exit_5", "exit_6", 
                                 "is_bike_la", "green_area", "poi_counts",
                                 "fictitious", "formOfWay", "density")]


# Min-Max scaling function
min_max_scaler <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Apply the scaler
highways_2020$length <- min_max_scaler(highways_2020$length)
highways_2020$slope <- min_max_scaler(highways_2020$slope)
highways_2020$green_area <- min_max_scaler(highways_2020$green_area)
highways_2020$poi_counts <- min_max_scaler(highways_2020$poi_counts)



summary(highways_2020$length)
summary(highways_2020$slope)
summary(highways_2020$exit_1)
summary(highways_2020$exit_2)
summary(highways_2020$exit_3)
summary(highways_2020$exit_4)
summary(highways_2020$exit_5)
summary(highways_2020$exit_6)
summary(highways_2020$is_bike_la)
summary(highways_2020$green_area)
summary(highways_2020$poi_counts)
summary(highways_2020$formOfWay)



crs <- st_crs(27700)
highways_2020 <- st_set_crs(highways_2020, crs)

#xmin <- 525000
#ymin <- 179000
#xmax <- 535000
#ymax <- 184000

xmin <- 527000
ymin <- 177000
xmax <- 537000
ymax <- 184000
bbox <- st_sf(geometry = st_sfc(st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax), c(xmin, ymin)))), crs = 27700))
clip_geometry <- st_as_sfc(bbox)
clipped_highways <- st_intersection(highways_2020, clip_geometry)


points <- st_coordinates(clipped_highways)
distances <- dist(points)
d <- 200  
W <- as.matrix(distances) < d
W <- 1 * W

clipped_highways$density_normalized <- (clipped_highways$density - min(clipped_highways$density)) / 
  (max(clipped_highways$density) - min(clipped_highways$density)) * 100
clipped_highways$density_normalized <- round(clipped_highways$density_normalized)


result <- S.CARbym(
  formula = density_normalized ~ length + slope + exit_1 + exit_2 + exit_3 + exit_4 + 
    exit_5 + exit_6 + is_bike_la + green_area + poi_counts + formOfWay,
  family = "poisson", 
  W = W, 
  data = clipped_highways, 
  burnin = 50000, 
  n.sample = 150000, 
  thin=100
)

summary(result)
summary(result$samples)
head(result$fitted.values)

head(clipped_highways$geometry)

clipped_highways$Fitted_Values <- result$fitted.values
print(head(clipped_highways))

new_df <- clipped_highways[, c("gml_id", "Fitted_Values")]
new_df$geometry <- NULL

head(new_df)
write.csv(new_df, "CAR_results_2020.csv", row.names = FALSE)





beta.samples <- result$samples$beta
library(coda)
plot(beta.samples[ ,2:4])


result$summary.results




hist(result$samples$beta, main = "Posterior distribution of beta", xlab = "Beta")

plot(clipped_highways$density_normalized, result$fitted.values, main = "Observed vs. fitted values", xlab = "Observed", ylab = "Fitted")

hist(residuals(result), main = "Histogram of residuals", xlab = "Residuals")

print(result$modelfit)

summary(result$samples$beta[, 1])
summary(result$samples$beta[, 2])
summary(result$samples$beta[, 3])
summary(result$samples$beta[, 4])
summary(result$samples$beta[, 5])
summary(result$samples$beta[, 6])
summary(result$samples$beta[, 7])
summary(result$samples$beta[, 8])
summary(result$samples$beta[, 9])
summary(result$samples$beta[, 10])
summary(result$samples$beta[, 11])
summary(result$samples$beta[, 12])
summary(result$samples$beta[, 13])







all_results <- data.frame(
  variables = c("length", "slope", "exit_1", "exit_2", "exit_3", "exit_4", "exit_5", "exit_6", "is_bike_la", "green_area", "poi_counts", "formOfWay"),
  Mean = c(-0.201011, -6.72202, -0.55471, -1.176097, -0.552112, -0.235596, 0.2957988, 0.455258, 0.550538, 0.5767822, -1.476610, 3.708223),
  SD = c(0.063711, 0.96972, 1.76743, 0.036188, 0.053399, 0.052384, 0.0291354, 0.061991, 0.084564, 0.0310894, 0.314195, 0.218004),
  Naive_SE = c(0.002015, 0.03067, 0.05589, 0.001144, 0.001689, 0.001657, 0.0009213, 0.001960, 0.002674, 0.0009831, 0.009936, 0.006894),
  Time_series_SE = c(0.004467, 0.04393, 0.06600, 0.002243, 0.004790, 0.004317, 0.0017880, 0.006132, 0.011886, 0.0017140, 0.016416, 0.012725),
  Quantile_2.5 = c(-0.3270, -8.550, -4.3585, -1.247, -0.6575, -0.3370, 0.2337, 0.3325, 0.3881, 0.5142, -2.0974, 3.298),
  Quantile_25 = c(-0.2419, -7.370, -1.6111, -1.200, -0.5855, -0.2696, 0.2768, 0.4129, 0.4897, 0.5559, -1.6879, 3.556),
  Quantile_50 = c(-0.2005, -6.734, -0.4127, -1.176, -0.5516, -0.2369, 0.2960, 0.4561, 0.5489, 0.5773, -1.4862, 3.705),
  Quantile_75 = c(-0.1587, -6.073, 0.6359, -1.153, -0.5173, -0.2033, 0.3166, 0.4978, 0.6054, 0.5980, -1.2577, 3.849),
  Quantile_97.5 = c(-0.0765, -4.876, 2.4084, -1.101, -0.4394, -0.1272, 0.3498, 0.5767, 0.7276, 0.6344, -0.8549, 4.147)
)

# 查看数据框
print(all_results)

