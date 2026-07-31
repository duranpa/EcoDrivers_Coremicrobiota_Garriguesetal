##Supplementary Figure 1b
#a) 100 sites across the South West of France
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
regions <- ne_states(country = "France", returnclass = "sf")
regions_mainland <- st_crop(regions, mainland_bbox)

# Crop France to mainland bbox
france_mainland <- st_crop(france, mainland_bbox)
regions_occ <- st_crop(regions, mainland_bbox)

# 4. Separate Occitanie and other regions
occitanie <- regions_mainland[regions_mainland$region %in% c("Occitanie"), ]
other_regions <- regions_mainland[!regions_mainland$region %in% c("Occitanie"), ]

# 5. Select your target region (replace "Bretagne" with your region name)
my_region <- regions_mainland %>% filter(name == "Occitanie")

sites$cluster_reduced <- factor(sites$cluster_reduced)

ggplot() +
  geom_sf(data = other_regions, fill = "#eeeeee", color = "black") +
  geom_sf(data = occitanie, fill = "#eeeeee", color = "black") +
  geom_sf(data = sites,
          aes(fill = cluster_reduced, alpha=0.5),
          shape = 21,          # supports fill + outline
          color = "black",     # stroke color
          size = 3,
          stroke = 0.3,        # stroke thickness
          alpha = 0.6) +
  scale_fill_manual(
    values = c(
      "incomplete_data" ="#d8d8d8",
      "1" = "#fb8072",
      "2" = "#6AA84F",
      "3" = "#47A5EC"
      )
  ) +
  theme_minimal() +
  labs(title = "Sites in Southwest France",
       fill = "Environmental clusters") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank()
  )

ggplot() +
  geom_sf(data = other_regions, fill = "#eeeeee", color = "black") +   # grey for others
  geom_sf(data = occitanie, fill = "#eeeeee", color = "black") +        # white for Occitanie
  ylim(42, 45)+xlim(-1, 4)+
  geom_sf(data = sites,
          aes(fill = cluster_reduced, alpha=0.5),
          shape = 21,          # supports fill + outline
          color = "black",     # stroke color
          size = 5,
          stroke = 0.3,        # stroke thickness
          alpha = 0.6) +
  scale_fill_manual(
    values = c(
      "incomplete_data" ="#d8d8d8",
      "1" = "#fb8072",
      "2" = "#6AA84F",
      "3" = "#47A5EC"
    )
  ) +              # sites colored by group
  theme_minimal() +
  labs(title = "Sites in Southwest France",
       fill = "Environmental clusters") +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank()
  )
