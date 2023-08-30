rm(list=ls())

library(sf)
library(CARBayes)



setwd("/Users/liangmingjun/Desktop/UCL/Individual Project")
# highways_centroid <- read.csv("highways_nkd_centroid.csv", row.names=NULL)
# file_path <- "Data/highways_nkd_centroid/highways_nkd_centroid.shp"
file_path <- "Data/final_highways/final_highways_2019.shp"
highways_2019 <- st_read(file_path)
names(highways_2019)
highways_2019 <- highways_2019[c("gml_id", "slope", "length", "exit_1", "exit_2", 
                                 "exit_3", "exit_4", "exit_5", "exit_6", 
                                 "is_bike_la", "green_area", "poi_counts",
                                 "fictitious", "formOfWay", "density")]


# Min-Max scaling function
min_max_scaler <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Apply the scaler
highways_2019$length <- min_max_scaler(highways_2019$length)
highways_2019$slope <- min_max_scaler(highways_2019$slope)
highways_2019$green_area <- min_max_scaler(highways_2019$green_area)
highways_2019$poi_counts <- min_max_scaler(highways_2019$poi_counts)



summary(highways_2019$length)
summary(highways_2019$slope)
summary(highways_2019$exit_1)
summary(highways_2019$exit_2)
summary(highways_2019$exit_3)
summary(highways_2019$exit_4)
summary(highways_2019$exit_5)
summary(highways_2019$exit_6)
summary(highways_2019$is_bike_la)
summary(highways_2019$green_area)
summary(highways_2019$poi_counts)
summary(highways_2019$formOfWay)



crs <- st_crs(27700)
highways_2019 <- st_set_crs(highways_2019, crs)

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
clipped_highways <- st_intersection(highways_2019, clip_geometry)


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

num_rows <- nrow(result$samples$beta)
print(num_rows)
result$samples$beta
summary(result$samples$beta[, 13])


summary(result)
summary(result$samples)
head(result$fitted.values)

head(clipped_highways$geometry)

library(sf)
library(leaflet)
clipped_highways$Fitted_Values <- result$fitted.values
print(head(clipped_highways))

# 提取 gml_id 和 Fitted_Values
new_df <- clipped_highways[, c("gml_id", "Fitted_Values")]
new_df$geometry <- NULL

# 显示前几行以确认数据
head(new_df)
# 导出数据到 CSV 文件
write.csv(new_df, "CAR_results_2019.csv", row.names = FALSE)





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
colnames(result$samples$beta)








all_results <- data.frame(
  variables = c("(Intercept)", "length", "slope", "exit_1", "exit_2", "exit_3", "exit_4", "exit_5", "exit_6", "is_bike_la", "green_area", "poi_counts", "formOfWay"),
  Mean = c(-0.174907, -5.92219, 0.55667, -1.114393, -0.558472, -0.151794, 0.225344, 0.465465, 0.398420, 0.608179, -0.31797, 2.850039, 0.0965431),
  SD = c(0.053954, 0.87917, 1.03374, 0.036108, 0.046036, 0.043231, 0.026626, 0.056514, 0.068790, 0.025867, 0.26848, 0.181149, 0.0083768),
  Naive_SE = c(0.001706, 0.02780, 0.03269, 0.001142, 0.001456, 0.001367, 0.000842, 0.001787, 0.002175, 0.000818, 0.00849, 0.005728, 0.0002649),
  Time_series_SE = c(0.003576, 0.03646, 0.03269, 0.002181, 0.003617, 0.002993, 0.001442, 0.004953, 0.007068, 0.001197, 0.01185, 0.011256, 0.0003223),
  Quantile_2.5 = c(-0.2811, -7.558, -1.5480, -1.187, -0.6557, -0.23828, 0.1722, 0.3455, 0.2575, 0.5605, -0.8460, 2.494, 0.07917),
  Quantile_25 = c(-0.2115, -6.538, -0.1219, -1.138, -0.5880, -0.18042, 0.2071, 0.4302, 0.3517, 0.5893, -0.4970, 2.729, 0.09118),
  Quantile_50 = c(-0.1739, -5.939, 0.5816, -1.114, -0.5587, -0.15079, 0.2253, 0.4666, 0.4022, 0.6083, -0.3211, 2.846, 0.09644),
  Quantile_75 = c(-0.1381, -5.339, 1.2643, -1.090, -0.5281, -0.12391, 0.2423, 0.5015, 0.4488, 0.6260, -0.1472, 2.974, 0.10227),
  Quantile_97.5 = c(-0.0723, -4.140, 2.4787, -1.043, -0.4724, -0.06848, 0.2768, 0.5746, 0.5213, 0.6574, 0.2207, 3.196, 0.11334)
)

# 查看数据框
print(all_results)

