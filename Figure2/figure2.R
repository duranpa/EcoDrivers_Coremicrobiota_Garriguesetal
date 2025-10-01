##figure 2
##figure2a
library(GUniFrac)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)

design<- read.table("design_noempty_bac.txt", header=TRUE, sep="\t")
asv_table<- read.table("asv_table_bac.txt", header=TRUE, sep="\t")
taxonomy<- read.table("taxonomy_mod_bac.txt", header=TRUE, sep="\t")

asv_raref <- Rarefy(t(asv_table), 2000)
asv_raref <- t(asv_raref$otu.tab.rff)
idx <- colnames(asv_table) %in% colnames(asv_raref)
asv_table <- asv_table[, idx]

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx <- taxonomy$Feature.ID %in% rownames(asv_table)
taxonomy<- taxonomy[idx,]

level <- which(colnames(taxonomy)=="core")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points_core <- points[!(points$Group.1 %in% c("non-core")),]

sd<- ddply(points_core, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
colors <- data.frame(group=c("Burkholderiales","Caulobacterales","Microtrichales",
                             "Propionibacteriales", "Rhizobiales","Solirubrobacterales","Sphingomonadales"),
                     color=c("#8bbfb1","#8ac095","#b9e0a3",
                             "#fffab4","#fcddb3","#b38e8e","#d9a28f"))


p<- ggplot(sd, aes(x = Site, y = mean, fill = Group.1))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  labs(y=paste("Average relative abundance(%)"), x=paste("Site"), fill="Core group")+
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.text = element_text(size=12), 
        legend.title = element_text(size=14), axis.title=element_text(size=14),)
p

#figure2b
design<- read.table("design_noempty_fun.txt", header=TRUE, sep="\t")
asv_table<- read.table("asv_table_fun.txt", header=TRUE, sep="\t")
taxonomy<- read.table("taxonomy_mod_fun.txt", header=TRUE, sep="\t")

asv_raref <- Rarefy(t(asv_table), 2000)
asv_raref <- t(asv_raref$otu.tab.rff)
idx <- colnames(asv_table) %in% colnames(asv_raref)
asv_table <- asv_table[, idx]

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

threshold <- .05
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx <- taxonomy$Feature.ID %in% rownames(asv_table)
taxonomy<- taxonomy[idx,]

level <- which(colnames(taxonomy)=="core")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points_core <- points[!(points$Group.1 %in% c("non-core")),]

sd<- ddply(points_core, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
colors <- data.frame(group=c("Chaetothyriales","Helotiales", "Hypocreales",
                             "Pleosporales","Xylariales"),
                     color=c("#737c64","#a7c0cf","#7c6473", 
                             "#d0a8bf","#bfd0a7"))

p<- ggplot(sd, aes(x = Site, y = mean, fill = Group.1))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  labs(y=paste("Average relative abundance(%)"), x=paste("Site"), fill="Core group")+
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.text = element_text(size=12), 
        legend.title = element_text(size=14), axis.title=element_text(size=14),)
p

#figure2c
number_sites<- read.table("number_ASVs_per_sitenew.txt", header=TRUE, sep="\t")

number_sites$core<- factor(number_sites$core, levels=c("Burkholderiales",
                                                         "Caulobacterales",
                                                         "Chaetothyriales",
                                                          "Helotiales",
                                                         "Hypocreales",
                                                         "Microtrichales",
                                                         "Pleosporales",
                                                         "Propionibacteriales",
                                                         "Rhizobiales",
                                                         "Solirubrobacterales",
                                                         "Sphingomonadales",
                                                         "Xylariales"))

colors <- data.frame(group=c("Burkholderiales","Caulobacterales","Chaetothyriales","Helotiales", "Hypocreales","Microtrichales","non-core",
                             "Pleosporales","Propionibacteriales", "Rhizobiales","Solirubrobacterales","Sphingomonadales",  
                             "Xylariales"),
                     color=c("#8bbfb1","#8ac095","#737c64","#a7c0cf","#7c6473", "#b9e0a3","#d9d9d9",
                             "#d0a8bf","#fffab4","#fcddb3","#b38e8e","#d9a28f","#bfd0a7"))


p<- ggplot(number_sites, aes(x = core, y = number_sites))+
  geom_violin()+
  geom_point(position = position_jitter(seed = 1, width = 0.3), color="black",alpha=0.2) +
  labs(y=paste("Number of sites (log10 transformed)"))+
  scale_y_log10()+
  geom_hline(yintercept=149,linetype=2)+
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
p+facet_grid(~Kingdom,scale="free_x", space="free")
