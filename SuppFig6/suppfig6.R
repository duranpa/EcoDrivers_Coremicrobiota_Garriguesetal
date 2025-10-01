##supplementary figure 6
library("factoextra")
library(FactoMineR)
env_factors<- read.table("environmental_factors.txt", header=TRUE, sep="\t")

climate_data <- scale(env_factors[, c( "Temp_high_Week1" ,                    "Temp_low_Week1",                      "Temp_avg_Week1",                      "Dewpoint_high_Week1" ,                "Dewpoint_low_Week1",                 
                                       "Dewpoint_avg_Week1",                  "Humidity_high_Week1" ,                "Humidity_low_Week1" ,                 "Humidity_avg_Week1" ,                 "Windgust_Week1",                     
                                       "Wind_high_Week1",                     "Wind_low_Week1",                      "Pressure_high_Week1",                 "Pressure_low_Week1" ,                 "Precipitation" ,                     
                                       "Temp_high_Week2",                     "Temp_low_Week2" ,                     "Temp_avg_Week2" ,                     "Dewpoint_high_Week2",                 "Dewpoint_low_Week2",                 
                                        "Dewpoint_avg_Week2",                  "Humidity_high_Week2",                 "Humidity_low_Week2" ,                 "Humidity_avg_Week2",                  "Windgust_Week2" ,                    
                                        "Wind_high_Week2" ,                    "Wind_low_Week2"  ,                    "Pressure_high_Week2" ,                "Pressure_low_Week2",                  "Precipitation.1",                    
                                        "Temp_high_Week3" ,                    "Temp_low_Week3"  ,                    "Temp_avg_Week3" ,                     "Dewpoint_high_Week3",                 "Dewpoint_low_Week3",                 
                                        "Dewpoint_avg_Week3" ,                 "Humidity_high_Week3",                 "Humidity_low_Week3",                  "Humidity_avg_Week3" ,                 "Windgust_Week3"  ,                   
                                        "Wind_high_Week3" ,                    "Wind_low_Week3"    ,                  "Pressure_high_Week3" ,                "Pressure_low_Week3" ,                 "Precipitation.2",                    
                                        "Temp_high_Week4" ,                    "Temp_low_Week4" ,                     "Temp_avg_Week4" ,                     "Dewpoint_high_Week4",                 "Dewpoint_low_Week4",                 
                                       "Dewpoint_avg_Week4",                  "Humidity_high_Week4",                 "Humidity_low_Week4",                  "Humidity_avg_Week4" ,                 "Windgust_Week4"  ,                   
                                       "Wind_high_Week4",                     "Wind_low_Week4" ,                     "Pressure_high_Week4",                 "Pressure_low_Week4",                  "Precipitation.3"       )])
soil_data <- scale(env_factors[, c("Total_nitrogen","Carbon_nitrogen_ratio","pH","Phosphore","Calcium","Magnesium","Sodium","Potassium","Iron","Aluminium","WHC","Organic_carbon","Soil_organic_matter","Manganese")])
plant_data <- scale(env_factors[, c("PCoA1_plant_communities",             "PCoA2_plant_communities" ,            "PCoA3_plant_communities" ,            "Shannon_diversity_plant_communities", "Plant_cover" )])

library(FactoMineR)

#plant_data <- plant_data[, apply(plant_data, 2, function(x) var(x, na.rm = TRUE) > 0)]

# Combine data

combined_data <- data.frame(climate_data, soil_data, plant_data)
group_sizes <- c(ncol(climate_data), ncol(soil_data), ncol(plant_data))

#combined_data <- combined_data[, apply(combined_data, 2, var) > 0]

mfa_res <- MFA(
  combined_data,
  group = group_sizes,
  type = c("s", "s", "s"),  # all scaled
  name.group = c("Climate", "Soil", "Plants")
)

site_coords <- mfa_res$ind$coord
dist_matrix <- dist(site_coords)
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc)

#to choose how many clusters
k2 <- kmeans(site_coords, centers = 3, nstart = 25)
str(k2)
fviz_cluster(k2, data = site_coords)

##elbow method to calculate the number of clusers
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(site_coords, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

#3 clusters
clusters <- cutree(hc, k = 4)

fviz_cluster(list(data = site_coords, cluster = clusters))
