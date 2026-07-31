##supplementary figure 6

library(GUniFrac)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)

##files are in folder of figure3
design_bac<- read.table("design_noempty_bac.txt", header=TRUE, sep="\t")
asv_table_bac<- read.table("asv_table_bac.txt", header=TRUE, sep="\t")
taxonomy_bac<- read.table("taxonomy_mod_bac.txt", header=TRUE, sep="\t")

design_fun<- read.table("design_noempty_fun.txt", header=TRUE, sep="\t")
asv_table_fun<- read.table("asv_table_fun.txt", header=TRUE, sep="\t")
taxonomy_fun<- read.table("taxonomy_mod_fun.txt", header=TRUE, sep="\t")

env_factors<- read.table("environmental_factors.txt", header=TRUE, sep="\t")

asv_raref_bac <- Rarefy(t(asv_table_bac), 2000)
asv_raref_bac <- t(asv_raref_bac$otu.tab.rff)
idx <- colnames(asv_table_bac) %in% colnames(asv_raref_bac)
asv_table_bac <- asv_table_bac[, idx]

asv_table_norm_bac <- apply(asv_table_bac, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm_bac * 100 > threshold) >= 1
asv_table_bac <- asv_table_bac[idx, ]
asv_table_norm_bac <- asv_table_norm_bac[idx, ]

idx <- taxonomy_bac$Feature.ID %in% rownames(asv_table_bac)
taxonomy_bac<- taxonomy_bac[idx,]

##RA of core ASVs
otu_long_bac <- melt(asv_table_norm_bac, id.vars = "Group.1", variable.name = "Sample")
points_bac <- cbind(otu_long_bac, design_bac[match(otu_long_bac$Var2, design_bac$SampleID), ])
points_bac_tax <- cbind(points_bac, taxonomy_bac[match(points_bac$Var1, taxonomy_bac$Feature.ID), ])

points_core_bac <- points_bac_tax[!(points_bac_tax$core %in% c("non-core")),]
points_core_bac2 <- cbind(points_core_bac, design_bac[match(points_core_bac$Var2, design_bac$SampleID), ])

points_core_bac3 <- points_core_bac2[, c("Alternate_sampleID", "Var1", "value")]
core_bac_table<- dcast(data = points_core_bac3,
                       formula = Alternate_sampleID ~ Var1,
                       value.var = "value")
core_bac_table2 <- merge(core_bac_table, design_bac[, c("Alternate_sampleID", "Site")],
                         by = "Alternate_sampleID",
                         all.x = TRUE)

asv_raref_fun <- Rarefy(t(asv_table_fun), 2000)
asv_raref_fun <- t(asv_raref_fun$otu.tab.rff)
idx <- colnames(asv_table_fun) %in% colnames(asv_raref_fun)
asv_table_fun <- asv_table_fun[, idx]

asv_table_norm_fun <- apply(asv_table_fun, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm_fun * 100 > threshold) >= 1
asv_table_fun <- asv_table_fun[idx, ]
asv_table_norm_fun <- asv_table_norm_fun[idx, ]

idx <- taxonomy_fun$Feature.ID %in% rownames(asv_table_fun)
taxonomy_fun<- taxonomy_fun[idx,]

otu_long_fun <- melt(asv_table_norm_fun, id.vars = "Group.1", variable.name = "Sample")
points_fun <- cbind(otu_long_fun, design_fun[match(otu_long_fun$Var2, design_fun$SampleID), ])
points_fun_tax <- cbind(points_fun, taxonomy_fun[match(points_fun$Var1, taxonomy_fun$Feature.ID), ])

points_core_fun <- points_fun_tax[!(points_fun_tax$core %in% c("non-core")),]

points_core_fun2 <- points_core_fun[, c("Alternate_sampleID", "Var1", "value")]
core_fun_table<- dcast(data = points_core_fun,
                       formula = Alternate_sampleID ~ Var1,
                       value.var = "value", fun.aggregate = mean)


#all together
points_core_all <- cbind(core_bac_table2, core_fun_table[match(core_bac_table2$Alternate_sampleID, core_fun_table$Alternate_sampleID), ])

points_core_all_env <- cbind(points_core_all, env_factors[match(points_core_all$Site, env_factors$Site), ])

points_core_all_env<- na.omit(points_core_all_env)
points_core_all_env_forcorr <- points_core_all_env[ , !(names(points_core_all_env) %in% c("Alternate_sampleID", "Site", "Alternate_sampleID", "Site")) ]

##since the table is too massive, I need to split it
set_env <- points_core_all_env_forcorr[, c("Total_nitrogen",	"Carbon_nitrogen_ratio",	"pH",	"Phosphore",	"Calcium",	"Magnesium",	"Sodium", 	"Potassium",	"Iron",	"Aluminium",	"WHC",	"Organic_carbon",	"Soil_organic_matter",	"Manganese",
                   "Temp_high_Week1",	"Temp_low_Week1",	"Temp_avg_Week1",	"Dewpoint_high_Week1",	"Dewpoint_low_Week1",	"Dewpoint_avg_Week1",	"Humidity_high_Week1",	"Humidity_low_Week1",	"Humidity_avg_Week1",	"Windgust_Week1", "Wind_high_Week1",	"Wind_low_Week1",	"Pressure_high_Week1",	"Pressure_low_Week1",	"Precipitation",
                   "Temp_high_Week2",	"Temp_low_Week2",	"Temp_avg_Week2",	"Dewpoint_high_Week2",	"Dewpoint_low_Week2",	"Dewpoint_avg_Week2",	"Humidity_high_Week2",	"Humidity_low_Week2",	"Humidity_avg_Week2",	"Windgust_Week2",	"Wind_high_Week2",	"Wind_low_Week2",	"Pressure_high_Week2",	"Pressure_low_Week2",	"Precipitation",
                   "Temp_high_Week3",	"Temp_low_Week3",	"Temp_avg_Week3",	"Dewpoint_high_Week3",	"Dewpoint_low_Week3",	"Dewpoint_avg_Week3",	"Humidity_high_Week3",	"Humidity_low_Week3",	"Humidity_avg_Week3",	"Windgust_Week3",	"Wind_high_Week3",	"Wind_low_Week3",	"Pressure_high_Week3",	"Pressure_low_Week3",	"Precipitation",	
                   "Temp_high_Week4",	"Temp_low_Week4",	"Temp_avg_Week4",	"Dewpoint_high_Week4",	"Dewpoint_low_Week4",	"Dewpoint_avg_Week4",	"Humidity_high_Week4",	"Humidity_low_Week4",	"Humidity_avg_Week4",	"Windgust_Week4",	"Wind_high_Week4",	"Wind_low_Week4",	"Pressure_high_Week4",	"Pressure_low_Week4",	"Precipitation",
                   "PCoA1_plant_communities",	"PCoA2_plant_communities",	"PCoA3_plant_communities",	"Shannon_diversity_plant_communities",	"Plant_cover")]
set_asvs <- points_core_all_env_forcorr[ , !(names(points_core_all_env_forcorr) %in% c("Total_nitrogen",	"Carbon_nitrogen_ratio",	"pH",	"Phosphore",	"Calcium",	"Magnesium",	"Sodium", 	"Potassium",	"Iron",	"Aluminium",	"WHC",	"Organic_carbon",	"Soil_organic_matter",	"Manganese",
                                                                                                          "Temp_high_Week1",	"Temp_low_Week1",	"Temp_avg_Week1",	"Dewpoint_high_Week1",	"Dewpoint_low_Week1",	"Dewpoint_avg_Week1",	"Humidity_high_Week1",	"Humidity_low_Week1",	"Humidity_avg_Week1",	"Windgust_Week1", "Wind_high_Week1",	"Wind_low_Week1",	"Pressure_high_Week1",	"Pressure_low_Week1",	"Precipitation",
                                                                                                          "Temp_high_Week2",	"Temp_low_Week2",	"Temp_avg_Week2",	"Dewpoint_high_Week2",	"Dewpoint_low_Week2",	"Dewpoint_avg_Week2",	"Humidity_high_Week2",	"Humidity_low_Week2",	"Humidity_avg_Week2",	"Windgust_Week2",	"Wind_high_Week2",	"Wind_low_Week2",	"Pressure_high_Week2",	"Pressure_low_Week2",	"Precipitation.1",
                                                                                                          "Temp_high_Week3",	"Temp_low_Week3",	"Temp_avg_Week3",	"Dewpoint_high_Week3",	"Dewpoint_low_Week3",	"Dewpoint_avg_Week3",	"Humidity_high_Week3",	"Humidity_low_Week3",	"Humidity_avg_Week3",	"Windgust_Week3",	"Wind_high_Week3",	"Wind_low_Week3",	"Pressure_high_Week3",	"Pressure_low_Week3",	"Precipitation.2",	
                                                                                                          "Temp_high_Week4",	"Temp_low_Week4",	"Temp_avg_Week4",	"Dewpoint_high_Week4",	"Dewpoint_low_Week4",	"Dewpoint_avg_Week4",	"Humidity_high_Week4",	"Humidity_low_Week4",	"Humidity_avg_Week4",	"Windgust_Week4",	"Wind_high_Week4",	"Wind_low_Week4",	"Pressure_high_Week4",	"Pressure_low_Week4",	"Precipitation.3",
                                                                                                          "PCoA1_plant_communities",	"PCoA2_plant_communities",	"PCoA3_plant_communities",	"Shannon_diversity_plant_communities",	"Plant_cover")) ]
##correlations
library(reshape2)
library(ggplot2)
library(vegan)
library(GUniFrac)
library(dplyr)
library(readr)
library(tidyverse)
library(ggsci)
library(vegan)
library(corrr)
library(corrplot)
library(Hmisc)
# Correlation matrix
cor_matrix <- cor(set_env, set_asvs, use = "pairwise.complete.obs", method = "spearman")

# P-value matrix
p.mat <- matrix(NA, nrow = ncol(set_env), ncol = ncol(set_asvs),
                dimnames = list(colnames(set_env), colnames(set_asvs)))

for (i in seq_len(ncol(set_env))) {
  for (j in seq_len(ncol(set_asvs))) {
    test <- cor.test(set_env[[i]], set_asvs[[j]], method = "spearman")
    p.mat[i, j] <- test$p.value
  }
}

p.vec <- as.vector(p.mat)
p.adj <- p.adjust(p.vec, method = "fdr")
p.adj.mat <- matrix(p.adj, nrow = nrow(p.mat), ncol = ncol(p.mat),
                    dimnames = dimnames(p.mat))


# Create the node and the edge
df1 <- melt(cor_matrix)
colnames(df1) <- c("node1", "node2", "edge")

df2 <- melt(p.adj.mat)
colnames(df2) <- c("node1", "node2", "sig")

df_net <- merge(df1, df2, by = c("node1", "node2"))
df_net<- na.omit(df_net)
df_net_sig <- df_net[df_net$sig < 0.05, ]
write.table(df_net_sig, "cor_core_asvs_RA_env.txt", sep="\t")
#modified table in excel to add metadata, tax assignment
cor_core_env<- read.table("cor_core_asvs_RA_env.txt", header=TRUE, sep="\t")

p<- ggplot(cor_core_env, aes(x = ASV, y = Factor, fill=edge))+
  geom_tile()+
  scale_fill_gradient2(
    low = "#950606", 
    mid = "#fffaa8", 
    high = "#069595", 
    midpoint = 0
  )+
  facet_grid(Category~core_group, scale="free", space="free")+
  labs(y=paste("Environmental factor"), x=paste("Core_group"), fill="Spearman's rho")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
p


