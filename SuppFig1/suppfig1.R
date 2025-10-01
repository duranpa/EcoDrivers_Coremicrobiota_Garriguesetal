##Supplementary Figure 1
#a) 121 sites across the South West of France
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(cowplot)

# France outline
france <- ne_countries(country = "France", scale = "large", returnclass = "sf")
# Load your sites
sites_df <- read.table("maps_clusters.txt", header=TRUE, sep="\t")  # must have columns: site, lon, lat
sites <- st_as_sf(sites_df, coords = c("longitude", "latitude"), crs = 4326)
sites$cluster_reduced<- as.factor(sites$cluster_reduced)
# Define mainland bounding box roughly
mainland_bbox <- st_bbox(c(xmin = -6, xmax = 10, ymin = 41, ymax = 52), crs = st_crs(france))
regions_mainland <- st_crop(regions, mainland_bbox)

# Crop France to mainland bbox
france_mainland <- st_crop(france, mainland_bbox)
regions_occ <- st_crop(regions, mainland_bbox)

# 4. Separate Occitanie and other regions
occitanie <- regions_mainland[regions_mainland$region %in% c("Occitanie"), ]
other_regions <- regions_mainland[!regions_mainland$region %in% c("Occitanie"), ]

# Regions (admin-1 boundaries)
regions <- ne_states(country = "France", returnclass = "sf")

# 4. Crop regions to mainland bbox for consistency
regions_mainland <- st_crop(regions, mainland_bbox)

# 5. Select your target region (replace "Bretagne" with your region name)
my_region <- regions_mainland %>% filter(name == "Occitanie")

ggplot() +
  geom_sf(data = other_regions, fill = "#eeeeee", color = "black") +   # grey for others
  geom_sf(data = occitanie, fill = "#eeeeee", color = "black") +        # white for Occitanie
  geom_sf(data = sites, aes(color = cluster_reduced, alpha=0.5), size = 3) +               # sites colored by group
  theme_minimal() +
  labs(title = "Sites in Southwest France",
       color = "Environmental clusters") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank()
  )

ggplot() +
  geom_sf(data = other_regions, fill = "#eeeeee", color = "black") +   # grey for others
  geom_sf(data = occitanie, fill = "#eeeeee", color = "black") +        # white for Occitanie
  ylim(42, 45)+xlim(-1, 4)+
  geom_sf(data = sites, aes(color = cluster_reduced, alpha=0.5), size = 3) +               # sites colored by group
  theme_minimal() +
  labs(title = "Sites in Southwest France",
       color = "Environmental clusters") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank()
  )

#b) climatic conditions
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

weather<- read.table("WU_metadata_long_subset.txt", header=TRUE, sep="\t")
sd<- ddply(weather, c("Factor", "Weeks", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

weather_scaled <- sd %>%
  group_by(Factor) %>%
  mutate(scaled_value = (mean - min(mean)) / (max(mean) - min(mean))) %>%
  ungroup()

# 1. Create a wide matrix for clustering
mat <- weather_scaled %>%
  unite("weeks_param", Weeks, Factor, sep = "_") %>%  # combine date + parameter into column id
  select(Site, weeks_param, scaled_value) %>%
  pivot_wider(names_from = weeks_param, values_from = scaled_value)

# 2. Make rownames = sites, and remove site column for distance calc
mat_sites <- as.data.frame(mat)
rownames(mat_sites) <- mat_sites$Site
mat_sites$Site <- NULL

# 3. Compute distance and perform clustering
d <- dist(mat_sites)                      # Euclidean distance
hc <- hclust(d, method = "ward.D2")       # hierarchical clustering

# 4. Get site order from clustering
site_order <- hc$labels[hc$order]

# 5. Reorder 'site' factor in original scaled dataset
weather_scaled$Site <- factor(weather_scaled$Site, levels = site_order)

# 6. Plot with reordered sites
ggplot(weather_scaled, aes(x = Weeks, y = Site, fill = scaled_value)) +
  geom_tile() +
  scale_fill_gradient(low = "grey", high = "black")+
  facet_grid(.~Factor) +
  theme_minimal() +
  theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1), axis.text.y = element_text(size=5)) +
  labs(x=paste("Weeks before harvesting"),
    fill = "Scaled value",
       title = "Normalized Climatic variables")

