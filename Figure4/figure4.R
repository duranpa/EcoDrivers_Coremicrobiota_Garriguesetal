##Figure3
#Figure3a
library(GUniFrac)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)

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

##RA of core orders
level <- which(colnames(taxonomy_bac)=="core")
asv_table_core_bac<- aggregate(asv_table_norm_bac, by=list(taxonomy_bac[, level]), FUN=sum)
otu_long_bac <- melt(asv_table_core_bac, id.vars = "Group.1", variable.name = "Sample")
points_bac <- cbind(otu_long_bac, design_bac[match(otu_long_bac$Sample, design_bac$SampleID), ])
points_core_bac <- points_bac[!(points_bac$Group.1 %in% c("non-core")),]
points_core_bac<- na.omit(points_core_bac)
points_core_bac2 <- points_core_bac[, c("Alternate_sampleID", "Group.1", "value")]
core_bac_table<- dcast(data = points_core_bac2,
                       formula = Alternate_sampleID ~ Group.1,
                       value.var = "value")
core_bac_table2 <- merge(core_bac_table, design_bac[, c("Alternate_sampleID", "Site")],
                   by = "Alternate_sampleID",
                   all.x = TRUE)



level <- which(colnames(taxonomy_fun)=="core")
asv_table_core_fun<- aggregate(asv_table_norm_fun, by=list(taxonomy_fun[, level]), FUN=sum)
otu_long_fun <- melt(asv_table_core_fun, id.vars = "Group.1", variable.name = "Sample")
points_fun <- cbind(otu_long_fun, design_fun[match(otu_long_fun$Sample, design_fun$SampleID), ])
points_core_fun <- points_fun[!(points_fun$Group.1 %in% c("non-core")),]
points_core_fun<- na.omit(points_core_fun)
points_core_fun2 <- points_core_fun[, c("Alternate_sampleID", "Group.1", "value")]
core_fun_table<- dcast(data = points_core_fun2,
                       formula = Alternate_sampleID ~ Group.1,
                       value.var = "value")
#all together
points_core_all <- cbind(core_bac_table2, core_fun_table[match(core_bac_table2$Alternate_sampleID, core_fun_table$Alternate_sampleID), ])

points_core_all_env <- cbind(points_core_all, env_factors[match(points_core_all$Site, env_factors$Site), ])

points_core_all_env<- na.omit(points_core_all_env)
points_core_all_env_forcorr <- points_core_all_env[ , !(names(points_core_all_env) %in% c("Alternate_sampleID", "Site", "Alternate_sampleID", "Site")) ]


#number of ASVs
otu_long_bac <- melt(asv_table_norm_bac, id.vars = "Group.1", variable.name = "Sample")
points_bac <- cbind(otu_long_bac, taxonomy_bac[match(otu_long_bac$Var1, taxonomy_bac$Feature.ID), ])

points_core_bac <- points_bac[!(points_bac$core %in% c("non-core")),]
points_core_bac2 <- cbind(points_core_bac, design_bac[match(points_core_bac$Var2, design_bac$SampleID), ])
df1_bac <- filter(points_core_bac2, value > 0)
sd_bac<- ddply(df1_bac, c("Site","core", "Alternate_sampleID"), summarise, N=length(value), sum=sum(value), sd=sd(value), se=sd/sqrt(N))
number_asvs_bac<- dcast(data = sd_bac,
                       formula = Alternate_sampleID+Site ~ core,
                       value.var = "N")

otu_long_fun <- melt(asv_table_norm_fun, id.vars = "Group.1", variable.name = "Sample")
points_fun <- cbind(otu_long_fun, taxonomy_fun[match(otu_long_fun$Var1, taxonomy_fun$Feature.ID), ])

points_core_fun <- points_fun[!(points_fun$core %in% c("non-core")),]
points_core_fun2 <- cbind(points_core_fun, design_fun[match(points_core_fun$Var2, design_fun$SampleID), ])
df1_fun <- filter(points_core_fun2, value > 0)
sd_fun<- ddply(df1_fun, c("core", "Alternate_sampleID"), summarise, N=length(value), sum=sum(value), sd=sd(value), se=sd/sqrt(N))
number_asvs_fun<- dcast(data = sd_fun,
                        formula = Alternate_sampleID ~ core,
                        value.var = "N")
#all together
number_asvs_all <- cbind(number_asvs_bac, number_asvs_fun[match(number_asvs_bac$Alternate_sampleID, number_asvs_fun$Alternate_sampleID), ])

number_asvs_all_env <- cbind(number_asvs_all, env_factors[match(number_asvs_all$Site, env_factors$Site), ])
number_asvs_all_env_forcorr <- number_asvs_all_env[ , !(names(number_asvs_all_env) %in% c("Alternate_sampleID", "Site", "Alternate_sampleID", "Site")) ]

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

M<- cor(number_asvs_all_env_forcorr, method='spearman', use="complete.obs") #or number of ASVs
testRes = cor.mtest(number_asvs_all_env_forcorr, conf.level = 0.95) #or number of ASVs
# Create the node and the edge
df1 <- melt(M)
colnames(df1) <- c("node1", "node2", "edge")

df2 <- melt(testRes$p)
colnames(df2) <- c("node1", "node2", "sig")

df_net <- merge(df1, df2, by = c("node1", "node2"))
df_net_sig <- df_net[df_net$sig < 0.05, ]
df_net_sig <- df_net_sig[df_net_sig$edge != 0, ]

#modified table in excel to contain categories and only one half of the correlation matrix, and both RA and number of ASVs
cor_core_env<- read.table("corr_core_number_asvs_envfactors.txt", header=TRUE, sep="\t")

p<- ggplot(cor_core_env, aes(x = Core_group_value, y = Factor, fill=edge))+
  geom_tile(color="black")+
  scale_fill_gradient2(
    low = "#950606", 
    mid = "#fffaa8", 
    high = "#069595", 
    midpoint = 0
  )+
  facet_grid(Category~Kingdom, scale="free", space="free")+
  labs(y=paste("Environmental factor"), x=paste("Core_group"), fill="Spearman's rho")+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1), 
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
p

#figure 3b

##RA of core ASVs
otu_long_bac <- melt(asv_table_norm_bac, id.vars = "Group.1", variable.name = "Sample")
points_bac <- cbind(otu_long_bac, design_bac[match(otu_long_bac$Var2, design_bac$SampleID), ])
points_bac_tax <- cbind(points_bac, taxonomy_bac[match(points_bac$Var1, taxonomy_bac$Feature.ID), ])

points_core_bac <- points_bac_tax[!(points_bac_tax$Correlation %in% c("not_correlated")),]


otu_long_fun <- melt(asv_table_norm_fun, id.vars = "Group.1", variable.name = "Sample")
points_fun <- cbind(otu_long_fun, design_fun[match(otu_long_fun$Var2, design_fun$SampleID), ])
points_fun_tax <- cbind(points_fun, taxonomy_fun[match(points_fun$Var1, taxonomy_fun$Feature.ID), ])

points_core_fun <- points_fun_tax[!(points_fun_tax$Correlation %in% c("not_correlated")),]

sum_bac<- ddply(points_core_bac, c("Correlation", "core", "Site","Correlation.1" ), summarise, N=length(value), sum=sum(value))
sd_bac<- ddply(sum_bac, c("Correlation", "core", "Correlation.1"), summarise, N=length(sum), mean=mean(sum), sd=sd(sum), se=sd/sqrt(N))

sum_fun<- ddply(points_core_fun, c("Correlation", "core", "Site",  "Correlation.1"), summarise, N=length(value), sum=sum(value))
sd_fun<- ddply(sum_fun, c("Correlation", "core",  "Correlation.1"), summarise, N=length(sum), mean=mean(sum), sd=sd(sum), se=sd/sqrt(N))

combined <- rbind(sd_bac, sd_fun)
##added kingdom
combined2<- read.table("RA_correlated_ASVs.txt", header=TRUE, sep="\t")
colors <- data.frame(group=c("Climatic conditions","Plant_diversity" , "Soil properties"),
                     color=c("#000004FF","#a32800" , "#FCFFA4FF"))

p<- ggplot(combined2, aes(x = Correlation, y = mean, fill=Category))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color))+
  labs(y=paste("Average relative abundance (%)"), x=paste("Spearman's correlation"))+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1), 
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
p+facet_grid(~Kingdom+core, space="free")

#or calculate proportion of ASVs
points_bac_tax <- cbind(points_bac, taxonomy_bac[match(points_bac$Var1, taxonomy_bac$Feature.ID), ])
points_bac_core <- points_bac_tax[!(points_bac_tax$core %in% c("non-core")),]

sum_bac<- ddply(points_bac_core, c("Correlation", "core", "Site","Correlation.1" ), summarise, N=length(value), sum=sum(value))

points_fun_tax <- cbind(points_fun, taxonomy_fun[match(points_fun$Var1, taxonomy_fun$Feature.ID), ])
points_fun_core <- points_fun_tax[!(points_fun_tax$core %in% c("non-core")),]

sum_fun<- ddply(points_fun_core, c("Correlation", "core", "Site","Correlation.1" ), summarise, N=length(value), sum=sum(value))

combined <- rbind(sum_bac, sum_fun)
#calculated proportion
combined3<- read.table("number_ASVs_correlated.txt", header=TRUE, sep="\t")
sd_combined<- ddply(combined3, c("Correlation", "core", "Category", "Kingdom"), summarise, N=length(proportion), mean=mean(proportion), sd=sd(proportion), se=sd/sqrt(N))

colors <- data.frame(group=c("Climatic conditions","Plant_diversity" , "Soil properties"),
                     color=c("#000004FF","#a32800" , "#FCFFA4FF"))

p<- ggplot(sd_combined, aes(x = Correlation, y = mean, fill=Category))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color))+
  labs(y=paste("Average proportion of correlated ASVs (%)"), x=paste("Spearman's correlation"))+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1), 
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
p+facet_grid(~Kingdom+core, space="free")

