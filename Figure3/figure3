##figure 2
##figure2a
library(GUniFrac)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)
library(grid)

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
points_core_bac <- points[!(points$Group.1 %in% c("non-core")),]

sd_bact<- ddply(points_core_bac, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
colors <- data.frame(group=c("Burkholderiales","Caulobacterales","Microtrichales",
                             "Propionibacteriales", "Rhizobiales","Solirubrobacterales","Sphingomonadales"),
                     color=c("#8bbfb1","#8ac095","#b9e0a3",
                             "#fffab4","#fcddb3","#b38e8e","#d9a28f"))


p<- ggplot(sd_bact, aes(x = Site, y = mean, fill = Group.1))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  labs(y=paste("Average relative abundance(%)"), x=paste("Site"), fill="Core group")+
  coord_flip() +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.text = element_text(size=12), 
        legend.title = element_text(size=14), axis.title=element_text(size=14),)
p

#figure2c
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
points_core_fun <- points[!(points$Group.1 %in% c("non-core")),]

sd_fun<- ddply(points_core_fun, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
colors <- data.frame(group=c("Chaetothyriales","Helotiales", "Hypocreales",
                             "Pleosporales","Xylariales"),
                     color=c("#737c64","#a7c0cf","#7c6473", 
                             "#d0a8bf","#bfd0a7"))

p<- ggplot(sd_fun, aes(x = Site, y = mean, fill = Group.1))+
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  labs(y=paste("Average relative abundance(%)"), x=paste("Site"), fill="Core group")+
  coord_flip() +
  theme(axis.text.y = element_blank(),
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


sd<- ddply(number_sites, c("core", "range", "Kingdom"), summarise, N=length(number_sites))
sd$range <- factor(sd$range ,
                     levels = c("80-90%",
                                "70-80%",
                                "60-70%",
                                "50-60%",
                                "40-50%",
                                "30-40%",
                                "20-30%",
                                "10-20%",
                                "5-10%",
                                "1-5%",
                                "0-1%"
                     ))

grad_cols <- colorRampPalette(c("#2986CC", "#e9f2f9"))(length(levels(sd$range)))
p<- ggplot(sd, aes(x = core, y = N, fill=range))+
  geom_bar(stat = "identity", colour="black") +
  scale_fill_viridis_d(option = "G")+
  labs(y=paste("Number of ASVs"))+
  scale_y_log10()+
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.text = element_text(size=12), 
  )
p+facet_grid(~Kingdom,scales="free_x", space="free")


p<- ggplot(sd, aes(x = core, y = range))+
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
