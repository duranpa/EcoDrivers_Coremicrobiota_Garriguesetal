##Figure1
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

##Figure1a
design<- read.table("design_noempty.txt", header=TRUE, sep="\t")
asv_table<- read.table("otu_table.txt", header=TRUE, sep="\t")
taxonomy<- read.table("taxonomy_mod.txt", header=TRUE, sep="\t")

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx <- taxonomy$Feature.ID %in% rownames(asv_table)
taxonomy<- taxonomy[idx,]

level <- which(colnames(taxonomy)=="Class")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points2 <- cbind(points,taxonomy[match(points$Group.1, taxonomy$Class), ])
points3 <- points2[points2$Algae_group %in% c("Green algae"), ]

sd<- ddply(points3, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
sd<- na.omit(sd)
p<- ggplot(sd, aes(x = reorder(Site, -mean), y = mean, fill=Group.1))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_brewer(palette = "Accent")+
  labs(y=paste("Relative abundance(%)"), x=paste("SW sites"),fill = "Algal class", size=16)+
  theme_bw()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##Figure1b
level <- which(colnames(taxonomy)=="Algae_group")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points3 <- points[points$Group.1 %in% c("Green algae"), ]
total_green_algae<- ddply(points3, c("Site"), summarise, N=length(value), mean=mean(value))

#bind with environmental factors
env_factors<- read.table("environmental_factors.txt", header=TRUE, sep="\t")
total_algae_env_factors <- cbind(points3,env_factors[match(points3$Site, env_factors$Site), ])
total_algae_env_factors<- na.omit(total_algae_env_factors)
write.table(total_algae_env_factors, "total_algae_env_factors_all.txt", sep="\t")
total_algae_env_factors<- read.table("total_algae_env_factors_all.txt", header=TRUE, sep="\t")

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

M<- cor(total_algae_env_factors, method='spearman', use="complete.obs")
testRes = cor.mtest(total_algae_env_factors, conf.level = 0.95)
# Create the node and the edge
df1 <- melt(M)
colnames(df1) <- c("node1", "node2", "edge")

df2 <- melt(testRes$p)
colnames(df2) <- c("node1", "node2", "sig")

df_net <- merge(df1, df2, by = c("node1", "node2"))
df_net_sig <- df_net[df_net$sig < 0.05, ]
df_net_sig <- df_net_sig[df_net_sig$edge != 0, ]

write.table(df_net_sig, "total_algae_env_factors_cor_sig.txt", sep="\t")
##modified table in excel for plotting

total_algae_env_factors_cor_sig<- read.table("total_algae_env_factors_cor_sig_plot.txt", header=TRUE, sep="\t")

colors <- data.frame(group=c("Climatic conditions","Plant_diversity" , "Soil properties"),
                     color=c("#000004FF","#a32800" , "#FCFFA4FF"))

p<- ggplot(total_algae_env_factors_cor_sig, aes(x = reorder(node2, -edge), y = edge, fill=Category))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color))+
  labs(y=paste("Spearman's rho"), x=paste("Environmental factor"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1), 
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
p

##Figure1c
level <- which(colnames(taxonomy)=="Class")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points2 <- cbind(points,taxonomy[match(points$Group.1, taxonomy$Class), ])
points3 <- points2[points2$Algae_group %in% c("Green algae"), ]
total_green_algae<- ddply(points3, c("Site", "Sample"), summarise, N=length(value), sum=sum(value))
avg_total_green_algae<- ddply(total_green_algae, c("Site"), summarise, N=length(sum), mean=mean(sum))
avg_algal_classes<- ddply(points3, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value))


algal_classes_env_factors <- cbind(avg_algal_classes,env_factors[match(avg_algal_classes$Site, env_factors$Site), ])
algal_classes_env_factors<- na.omit(algal_classes_env_factors)
write.table(algal_classes_env_factors, "algal_classes_env_factors.txt", sep="\t")
algal_classes_env_factors<- read.table("algal_classes_env_factors.txt", header=TRUE, sep="\t")

M<- cor(algal_classes_env_factors, method='spearman', use="complete.obs")
testRes = cor.mtest(algal_classes_env_factors, conf.level = 0.95)
# Create the node and the edge
df1 <- melt(M)
colnames(df1) <- c("node1", "node2", "edge")

df2 <- melt(testRes$p)
colnames(df2) <- c("node1", "node2", "sig")

df_net <- merge(df1, df2, by = c("node1", "node2"))
df_net_sig <- df_net[df_net$sig < 0.05, ]
df_net_sig <- df_net_sig[df_net_sig$edge != 0, ]
write.table(df_net_sig, "algal_classes_env_factors_cor_sig.txt", sep="\t")
##modified table in excel for plotting

algal_classes_env_factors_cor_sig<- read.table("algal_classes_env_factors_cor_sig_plot.txt", header=TRUE, sep="\t")

p<- ggplot(algal_classes_env_factors_cor_sig, aes(x = Factor, y = Algal_class, fill=edge))+
  geom_tile()+
  scale_fill_gradient2(
    low = "#950606", 
    mid = "#fffaa8", 
    high = "#069595", 
    midpoint = 0
  )+
  facet_grid(.~Category, scale="free_x", space="free_x")+
  labs(y=paste("Algal class"), x=paste("Environmental factor"), fill="Spearman's rho")+
  theme_bw()+
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1), 
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))
p
