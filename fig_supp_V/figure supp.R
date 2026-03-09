
library(fmsb)
library(dplyr)
library(tidyr)
library(vegan)
library(tidyverse)

## load data
cluster <- read.table("20250818_clusters_env_factors.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)
cluster <- cluster[,-1]

taxonomy_bac <- read.table("taxonomy_mod_bac.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)
taxonomy_fun <- read.table("taxonomy_mod_fun.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)

taxonomy_bac  <- taxonomy_bac  %>% rename(Kingdom = Domain)
cols_to_keep <- c("Feature.ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy_fun <- taxonomy_fun[, cols_to_keep]
taxonomy_bac <- taxonomy_bac[, cols_to_keep]

taxonomy_all <- rbind(taxonomy_bac , taxonomy_fun)



## ASVs
asv_table_bac <- read.table("asv_table_bac.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)
asv_table_fun <- read.table("asv_table_fun.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)

asv_raref_bac <- rrarefy(t(asv_table_bac), 2000)
asv_raref_bac <- t(asv_raref_bac)
idx <- colnames(asv_table_bac) %in% colnames(asv_raref_bac)
asv_table_bac <- asv_table_bac[, idx]

asv_table_norm_bac <- apply(asv_table_bac, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm_bac * 100 > threshold) >= 1
asv_table_bac <- asv_table_bac[idx, ]
asv_table_norm_bac <- asv_table_norm_bac[idx, ]

idx <- taxonomy_bac$Feature.ID %in% rownames(asv_table_bac)
taxonomy_bac<- taxonomy_bac[idx,]

asv_raref_fun <- rrarefy(t(asv_table_fun), 2000)
asv_raref_fun <- t(asv_raref_fun)
idx <- colnames(asv_table_fun) %in% colnames(asv_raref_fun)
asv_table_fun <- asv_table_fun[, idx]

asv_table_norm_fun <- apply(asv_table_fun, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm_fun * 100 > threshold) >= 1
asv_table_fun <- asv_table_fun[idx, ]
asv_table_norm_fun <- asv_table_norm_fun[idx, ]

idx <- taxonomy_fun$Feature.ID %in% rownames(asv_table_fun)
taxonomy_fun<- taxonomy_fun[idx,]


# modifier, ajouter les clusters et merge les design file bac et fun
design_fun<- read.table("design_noempty_fun.txt", header=TRUE, sep="\t")
design_bac<- read.table("design_noempty_bac.txt", header=TRUE, sep="\t")

design_bac_new <- merge(design_bac, cluster[, c("Site.1", "Cluster")],
                        by.x = "Site", 
                        by.y = "Site.1",
                        all.x = TRUE)

design_fun_new <- merge(design_fun, cluster[, c("Site.1", "Cluster")],
                        by.x = "Site", 
                        by.y = "Site.1",
                        all.x = TRUE)


## merged les fichiers ASV et NA --> 0
ASVs <- bind_rows(as.data.frame(asv_table_norm_bac),as.data.frame(asv_table_norm_fun))
ASVs[is.na(ASVs)] <- 0


## separe les infos de desing de chaque cluster 
design_bac_new <- design_bac_new[,-6]

design_all <- rbind(design_bac_new , design_fun_new)
design_all <- na.omit(design_all)

design_cluster1 <- design_all[design_all$Cluster == "1", ]
design_cluster2 <- design_all[design_all$Cluster == "2", ]
design_cluster3 <- design_all[design_all$Cluster == "3", ]


## garder que les ASVs present dans cluster 1/2/3 séparement pour faire 3 reseaux a part

## Extraire la liste des SampleID à garder
samples_to_keep1 <- design_cluster1$SampleID
samples_to_keep2 <- design_cluster2$SampleID
samples_to_keep3 <- design_cluster3$SampleID

## Garder uniquement les colonnes correspondantes
ASVs_cluster1 <- ASVs[, colnames(ASVs) %in% samples_to_keep1]
ASVs_cluster2 <- ASVs[, colnames(ASVs) %in% samples_to_keep2]
ASVs_cluster3 <- ASVs[, colnames(ASVs) %in% samples_to_keep3]

ASVs_cluster1 <- ASVs_cluster1[rowSums(ASVs_cluster1) != 0, ]
ASVs_cluster1 <- ASVs_cluster1 %>% rownames_to_column(var = "Feature.ID")

ASVs_cluster2 <- ASVs_cluster2[rowSums(ASVs_cluster2) != 0, ]
ASVs_cluster2 <- ASVs_cluster2 %>% rownames_to_column(var = "Feature.ID")

ASVs_cluster3 <- ASVs_cluster3[rowSums(ASVs_cluster3) != 0, ]
ASVs_cluster3 <- ASVs_cluster3 %>% rownames_to_column(var = "Feature.ID")

taxonomy_all <- rbind(taxonomy_bac , taxonomy_fun)




#### NETWORK ANALYSIS ####

node_C1_F <- read_csv("node_C1_F.csv")
node_C2_F <- read_csv("node_C2_F.csv")
node_C3_F <- read_csv("node_C3_F.csv")


group <- c("Burkholderiales", "Caulobacterales", "Frankiales", "Microtrichales",
           "Propionibacteriales", "Helotiales", "Hypocreales", "Pleosporales",
           "Rhizobiales", "Sphingomonadales", "Xanthomonadales", "Chaetothyriales",
           "Solirubrobacterales", "Chloroflexales", "Gemmatimonadales", "Micrococcales",
           "Pseudonocardiales", "Xylariales")

color <- c("#B5EAD7", "#FFDAC1", "#FF9AA2", "#C7CEEA",
           "#FFB7B2", "#F1F0C0", "#B0E0A8", "#F6C6EA",
           "#C4F0C5", "#E6D0DE", "#FFD3B6", "#E0BBE4", "#D1F5D3",
           "#FFE0AC", "#E2F0CB", "#D3E5FF", "#FDE2E4", "#e3daff" )

node_C1_F$Core_status <- ifelse(node_C1_F$Order %in% group,
                                "Core",
                                "Non_Core")
node_C2_F$Core_status <- ifelse(node_C2_F$Order %in% group,
                                "Core",
                                "Non_Core")
node_C3_F$Core_status <- ifelse(node_C3_F$Order %in% group,
                                "Core",
                                "Non_Core")


##### plot pour core vs non core
node_scaled <- node_C3_F %>%
       mutate(across(
             c(BetweennessCentrality, Degree, ClosenessCentrality),
             ~ (. - min(., na.rm = TRUE)) /
                   (max(., na.rm = TRUE) - min(., na.rm = TRUE)),
             .names = "{.col}_scaled"
         ))

summary_scaled <- node_scaled %>%
       group_by(Core_status) %>%
       summarise(
            Betweenness = mean(BetweennessCentrality_scaled, na.rm = TRUE),
            Degree      = mean(Degree_scaled, na.rm = TRUE),
            Closeness   = mean(ClosenessCentrality_scaled, na.rm = TRUE)
         )


radar_long <- summary_scaled %>%
  pivot_longer(
    cols = -Core_status,
    names_to = "Metric",
    values_to = "Value"
  )

radar_long$Metric <- factor(
  radar_long$Metric,
  levels = c("Betweenness", "Degree", "Closeness")
)


# Transform in wide
radar_wide <- radar_long %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  column_to_rownames("Core_status")

# Add max & min for fmsb
radar_data <- rbind(
  rep(1, ncol(radar_wide)),  # max
  rep(0, ncol(radar_wide)),  # min
  radar_wide
)

# Colors
colors <- c("red", "blue")  # Core, Non_Core

# Plot
radarchart(radar_data, 
           pcol = colors, 
           pfcol = scales::alpha(colors, 0.3), 
           plwd = 2,
           title = "Radar chart Core vs Non core (cluster 3)",
           cglcol = "grey", 
           cglty = 1, 
           axislabcol = "grey", 
           vlcex = 1.1)
legend("topright", legend = rownames(radar_wide), col = colors, lty = 1, lwd = 2)



## by core members


library(dplyr)
library(tidyr)
library(fmsb)
library(scales)
library(tibble)

group <- c("Burkholderiales", "Caulobacterales", "Frankiales", "Microtrichales",
           "Propionibacteriales", "Helotiales", "Hypocreales", "Pleosporales",
           "Rhizobiales", "Sphingomonadales", "Xanthomonadales", "Chaetothyriales",
           "Solirubrobacterales", "Chloroflexales", "Gemmatimonadales", "Micrococcales",
           "Pseudonocardiales", "Xylariales")

color <- c("#B5EAD7", "#FFDAC1", "#FF9AA2", "#C7CEEA",
           "#FFB7B2", "#F1F0C0", "#B0E0A8", "#F6C6EA",
           "#C4F0C5", "#E6D0DE", "#FFD3B6", "#E0BBE4", "#D1F5D3",
           "#FFE0AC", "#E2F0CB", "#D3E5FF", "#FDE2E4", "#e3daff")

selected_orders <- c("Burkholderiales", "Caulobacterales", "Microtrichales",
                     "Propionibacteriales", "Rhizobiales", "Solirubrobacterales",
                     "Sphingomonadales", "Chaetothyriales", "Helotiales",
                     "Hypocreales", "Pleosporales", "Xylariales")

color_map <- setNames(color, group)

# Normalisation
node_scaled <- node_C1_F %>%
  mutate(
    Betweenness = rescale(BetweennessCentrality),
    Degree = rescale(Degree),
    Closeness = rescale(ClosenessCentrality)
  )

# Mean by Order 
core_summary <- node_scaled %>%
  filter(Core_status == "Core", Order %in% selected_orders) %>%
  group_by(Order) %>%
  summarise(
    Betweenness = mean(Betweenness),
    Degree = mean(Degree),
    Closeness = mean(Closeness),
    .groups = "drop"
  )

# Mean Non_Core
noncore_summary <- node_scaled %>%
  filter(Core_status == "Non_Core") %>%
  summarise(
    Betweenness = mean(Betweenness),
    Degree = mean(Degree),
    Closeness = mean(Closeness)
  ) %>%
  mutate(Order = "Non_Core")

# Transform in wide
radar_wide <- bind_rows(core_summary, noncore_summary) %>%
  column_to_rownames("Order")

# Add max & min for fmsb
radar_data <- rbind(
  rep(1, ncol(radar_wide)),
  rep(0, ncol(radar_wide)),
  radar_wide
)


plot_colors <- c(color_map[rownames(radar_wide)[-nrow(radar_wide)]], "black")

# Plot
par(mar = c(1,1,3,1))
radarchart(radar_data,
           pcol = plot_colors,
           plwd = 2,
           plty = 1,
           title = "Radar chart Core vs Non core for cluster 1",
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           vlcex = 1.1)

legend("topright",
       legend = rownames(radar_wide),
       col = plot_colors,
       lty = 1,
       lwd = 2,
       cex = 0.6,
       ncol = 2,
       bty = "n")
