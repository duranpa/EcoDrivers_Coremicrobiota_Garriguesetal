
## Abundance - Occupancy analysis based on Shade & Stopnisek (2019) ##

library(ggplot2)
library(dplyr)
library(tidyverse)
library(GUniFrac)
library(ggrepel)
library(readxl)
library(vegan)

## Load data 
design_bac <- read_excel("data_occupancy/design_all_bac.xlsx")

asv_table_bac <- read.table("data_occupancy/occ_all_bac2.txt", header=T, sep="\t", check.names = F, row.names = NULL)

# Photosynthetic - microbiota 
asv_table_bac <- asv_table_bac %>% column_to_rownames(var = "Group.1")
asv_table_bac <- t(asv_table_bac)
asv_table_bac <- as.data.frame(asv_table_bac)

nReads=1000 

# occupancy
otu <- asv_table_bac
otu_PA <- 1 * ((otu > 0) == 1)                                                  # Présence-absence data
otu_occ <- rowSums(otu_PA) / ncol(otu_PA)  # Occupancy calculation
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)         # Mean relative abundance

occ_abun_bac <- data.frame(otu_occ = otu_occ, otu_rel = otu_rel) %>%
  rownames_to_column('Group.1')


# Ranking OTUs based on their occupancy
# For calculating ranking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
design_bac_unique <- design_bac %>% distinct(s_h, .keep_all = TRUE)


PresenceSum <- data.frame(Group.1 = as.factor(row.names(otu)), otu) %>%             # Modify Photsynthetic - Microbiota 
  gather(s_h, abun, -Group.1) %>%                                              # - design_bac (bac) + s_h
  left_join(design_bac_unique , by = 's_h') %>%                                 # - design_fungi (fungi) + s_h
  group_by(Group.1, s_h) %>%
  summarise(freq_x = sum(abun > 0) / length(abun),                              
            core_x = ifelse(freq_x == 1, 1, 0)) %>%                             
  group_by(Group.1) %>%                                                             
  summarise(sumF = sum(freq_x),
            sumG = sum(core_x),
            nS = length(s_h),
            Index = (sumF + sumG) / nS)  

otu_ranked <- occ_abun_bac %>%
  left_join(PresenceSum, by = 'Group.1') %>%
  transmute(Group.1 = Group.1,
            rank = Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

# calculating BC dissimilarity based on the 1st ranked OTU
otu_start <- otu_ranked$Group.1[1]  
start_matrix <- as.matrix(otu[otu_start,])
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
df_s <- data.frame(x_names, x)
names(df_s)[2] <- 1
BCaddition <- rbind(BCaddition, df_s)
# calculating BC dissimilarity based on addition of ranked OTUs from 2nd to 500th. 
# Can be set to the entire length of OTUs in the dataset, however it might take 
# some time if more than 5000 OTUs are included.

for (i in 2:length(otu_ranked$Group.1)) {
  otu_add <- otu_ranked$Group.1[i]
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
  df_a <- data.frame(x_names, x)
  names(df_a)[2] <- i
  BCaddition <- left_join(BCaddition, df_a, by = 'x_names')
}

# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>%
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC = mean(BC, na.rm = TRUE)) %>%                                # Mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC = MeanBC / max(MeanBC, na.rm = TRUE))                     

if ("rank.x" %in% names(BC_ranked) & "rank.y" %in% names(BC_ranked)) {
  BC_ranked <- BC_ranked %>%
    select(-rank.y) %>%   
    rename(rank = rank.x)
}

# Modify for each taxnomic level:
# 182.x & 182.y for Photosynthetic - Microbiota (fungi)
# 543.x & 543.y for Photosynthetic - Microbiota (bac)

if ("543.x" %in% BC_ranked$rank | "543.y" %in% BC_ranked$rank) {
  BC_ranked$rank <- gsub("\\.x|\\.y", "", BC_ranked$rank)  
  BC_ranked$rank <- as.numeric(BC_ranked$rank)  
}

Increase <- BC_ranked$MeanBC[-1] / BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0, Increase), rank = 1:(length(Increase) + 1))
increaseDF$rank<- as.factor(increaseDF$rank)
BC_ranked$rank<- as.factor(BC_ranked$rank)

BC_ranked <- left_join(BC_ranked, increaseDF, by = "rank")
BC_ranked <- BC_ranked[-nrow(BC_ranked),]


#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])))

# Create a column defining "core" OTUs without Unassigned & 0
occ_abun_bac <- occ_abun_bac[occ_abun_bac$Group.1 != "Unassigned", ]
occ_abun_bac <- occ_abun_bac[occ_abun_bac$Group.1 != "0", ]
occ_abun_bac$fill <- 'no'
occ_abun_bac$fill[occ_abun_bac$Group.1 %in% otu_ranked$Group.1[1:elbow]] <- 'core'


group <- c("Burkholderiales","Caulobacterales","Chaetothyriales","Helotiales", "Hypocreales","Microtrichales","non-core",
           "Pleosporales","Propionibacteriales", "Rhizobiales","Solirubrobacterales","Sphingomonadales",  
           "Xylariales")
color <- c("#8bbfb1","#8ac095","#737c64","#a7c0cf","#7c6473", "#b9e0a3","#d9d9d9",
           "#d0a8bf","#fffab4","#fcddb3","#b38e8e","#d9a28f","#bfd0a7" )


colors_df <- data.frame(
  group = group,
  color = color,
  stringsAsFactors = FALSE)

elbow_points <- occ_abun_bac[occ_abun_bac$Group.1 %in% otu_ranked$Group.1[1:elbow], ]

occ_abun_bac <- occ_abun_bac[, !duplicated(colnames(occ_abun_bac))]
occ_abun_bac$Group.1 <- as.character(occ_abun_bac$Group.1)
colors_df$group <- as.character(colors_df$group)

occ_abun <- merge(
  occ_abun_bac[, c("Group.1", "otu_occ", "otu_rel", "fill")],
  colors_df,
  by.x = "Group.1", by.y = "group",
  all.x = TRUE)


# Plot
p1 <- ggplot() +
  geom_point(data = subset(occ_abun, fill == "no"), aes(x = log10(otu_rel), y = otu_occ), pch = 21, fill = "grey", alpha = 0.5, size = 2.5) +
  geom_point(data = subset(occ_abun, fill == "core" & !is.na(color)), aes(x = log10(otu_rel), y = otu_occ, fill = color), pch = 21, size = 4, color = "black", show.legend = FALSE) +
  geom_text_repel(data = subset(occ_abun, fill == "core" & !is.na(color)), aes(x = log10(otu_rel), y = otu_occ, label = Group.1),
                  size = 5, color = "black", box.padding = 1.2, point.padding = 0.8, force = 4, force_pull = 0.2, max.iter = 10000, max.overlaps = Inf, segment.color = "black", segment.size = 0.3, show.legend = FALSE) +
  geom_vline(xintercept = log10(min(elbow_points$otu_rel)), color = "red", linetype = "dashed") +
  scale_fill_identity() +
  labs(x = "Relative abundance (log10)", y = "Occupancy") +
  coord_cartesian(clip = "off") +
  expand_limits(y = 1.1) +  
  theme_classic() +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), plot.margin = margin(t = 20, r = 80, b = 10, l = 10))
p1

##fungi

## Load data 
design_fungi <- read_excel("data_occupancy/design_all_fungi.xlsx")

asv_table_fun <- read.table("data_occupancy/occ_all_fg.txt", header=T, sep="\t", check.names = F)


# Photosynthetic - microbiota 
asv_table_fun <- asv_table_fun %>% column_to_rownames(var = "Group.1")
asv_table_fun <- t(asv_table_fun)
asv_table_fun <- as.data.frame(asv_table_fun)

nReads=1000 

# occupancy
otu <- asv_table_fun
otu_PA <- 1 * ((otu > 0) == 1)                                                  # Présence-absence data
otu_occ <- rowSums(otu_PA) / ncol(otu_PA)  # Occupancy calculation
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)         # Mean relative abundance

occ_abun_fun <- data.frame(otu_occ = otu_occ, otu_rel = otu_rel) %>%
  rownames_to_column('Group.1')


# Ranking OTUs based on their occupancy
# For calculating ranking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
design_fun_unique <- design_fungi %>% distinct(s_h, .keep_all = TRUE)


PresenceSum <- data.frame(Group.1 = as.factor(row.names(otu)), otu) %>%             # Modify Photsynthetic - Microbiota 
  gather(s_h, abun, -Group.1) %>%                                              # - design_bac (bac) + s_h
  left_join(design_fun_unique , by = 's_h') %>%                                 # - design_fungi (fungi) + s_h
  group_by(Group.1, s_h) %>%
  summarise(freq_x = sum(abun > 0) / length(abun),                              
            core_x = ifelse(freq_x == 1, 1, 0)) %>%                             
  group_by(Group.1) %>%                                                             
  summarise(sumF = sum(freq_x),
            sumG = sum(core_x),
            nS = length(s_h),
            Index = (sumF + sumG) / nS)  

otu_ranked <- occ_abun_fun %>%
  left_join(PresenceSum, by = 'Group.1') %>%
  transmute(Group.1 = Group.1,
            rank = Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

# calculating BC dissimilarity based on the 1st ranked OTU
otu_start <- otu_ranked$Group.1[1]  
start_matrix <- as.matrix(otu[otu_start,])
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
df_s <- data.frame(x_names, x)
names(df_s)[2] <- 1
BCaddition <- rbind(BCaddition, df_s)
# calculating BC dissimilarity based on addition of ranked OTUs from 2nd to 500th. 
# Can be set to the entire length of OTUs in the dataset, however it might take 
# some time if more than 5000 OTUs are included.

for (i in 2:length(otu_ranked$Group.1)) {
  otu_add <- otu_ranked$Group.1[i]
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
  df_a <- data.frame(x_names, x)
  names(df_a)[2] <- i
  BCaddition <- left_join(BCaddition, df_a, by = 'x_names')
}

# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>%
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC = mean(BC, na.rm = TRUE)) %>%                                # Mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC = MeanBC / max(MeanBC, na.rm = TRUE))                     

if ("rank.x" %in% names(BC_ranked) & "rank.y" %in% names(BC_ranked)) {
  BC_ranked <- BC_ranked %>%
    select(-rank.y) %>%   
    rename(rank = rank.x)
}

# Modify for each taxnomic level:
# 182.x & 182.y for Photosynthetic - Microbiota (fungi)
# 543.x & 543.y for Photosynthetic - Microbiota (bac)

if ("182.x" %in% BC_ranked$rank | "182.y" %in% BC_ranked$rank) {
  BC_ranked$rank <- gsub("\\.x|\\.y", "", BC_ranked$rank)  
  BC_ranked$rank <- as.numeric(BC_ranked$rank)  
}

Increase <- BC_ranked$MeanBC[-1] / BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0, Increase), rank = 1:(length(Increase) + 1))
increaseDF$rank<- as.factor(increaseDF$rank)
BC_ranked$rank<- as.factor(BC_ranked$rank)

BC_ranked <- left_join(BC_ranked, increaseDF, by = "rank")
BC_ranked <- BC_ranked[-nrow(BC_ranked),]


#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])))

# Create a column defining "core" OTUs without Unassigned & 0
occ_abun_fun <- occ_abun_fun[occ_abun_fun$Group.1 != "Unassigned", ]
occ_abun_fun <- occ_abun_fun[occ_abun_fun$Group.1 != "0", ]
occ_abun_fun$fill <- 'no'
occ_abun_fun$fill[occ_abun_fun$Group.1 %in% otu_ranked$Group.1[1:elbow]] <- 'core'


group <- c("Burkholderiales","Caulobacterales","Chaetothyriales","Helotiales", "Hypocreales","Microtrichales","non-core",
           "Pleosporales","Propionibacteriales", "Rhizobiales","Solirubrobacterales","Sphingomonadales",  
           "Xylariales")
color <- c("#8bbfb1","#8ac095","#737c64","#a7c0cf","#7c6473", "#b9e0a3","#d9d9d9",
           "#d0a8bf","#fffab4","#fcddb3","#b38e8e","#d9a28f","#bfd0a7" )


colors_df <- data.frame(
  group = group,
  color = color,
  stringsAsFactors = FALSE)

elbow_points <- occ_abun_fun[occ_abun_fun$Group.1 %in% otu_ranked$Group.1[1:elbow], ]

occ_abun_fun <- occ_abun_fun[, !duplicated(colnames(occ_abun_fun))]
occ_abun_fun$Group.1 <- as.character(occ_abun_fun$Group.1)
colors_df$group <- as.character(colors_df$group)

occ_abun <- merge(
  occ_abun_fun[, c("Group.1", "otu_occ", "otu_rel", "fill")],
  colors_df,
  by.x = "Group.1", by.y = "group",
  all.x = TRUE)


# Plot
p1 <- ggplot() +
  geom_point(data = subset(occ_abun, fill == "no"), aes(x = log10(otu_rel), y = otu_occ), pch = 21, fill = "grey", alpha = 0.5, size = 2.5) +
  geom_point(data = subset(occ_abun, fill == "core" & !is.na(color)), aes(x = log10(otu_rel), y = otu_occ, fill = color), pch = 21, size = 4, color = "black", show.legend = FALSE) +
  geom_text_repel(data = subset(occ_abun, fill == "core" & !is.na(color)), aes(x = log10(otu_rel), y = otu_occ, label = Group.1),
                  size = 5, color = "black", box.padding = 1.2, point.padding = 0.8, force = 4, force_pull = 0.2, max.iter = 10000, max.overlaps = Inf, segment.color = "black", segment.size = 0.3, show.legend = FALSE) +
  geom_vline(xintercept = log10(min(elbow_points$otu_rel)), color = "red", linetype = "dashed") +
  scale_fill_identity() +
  labs(x = "Relative abundance (log10)", y = "Occupancy") +
  coord_cartesian(clip = "off") +
  expand_limits(y = 1.1) +  
  theme_classic() +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), plot.margin = margin(t = 20, r = 80, b = 10, l = 10))
p1

