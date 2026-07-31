##supplementary figure 3
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

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

level <- which(colnames(taxonomy)=="Phylum")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])

#average relative abundance across all sites
sd<- ddply(points, c("Group.1"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
sd1 <- cbind(sd, taxonomy[match(sd$Group.1, taxonomy$Phylum), ])

colors <- data.frame(group=c("Green algae", "Land plants", "Other eukaryotes"),
                     color=c("#32a852", "#38761d", "#808080"))

p<- ggplot(sd1, aes(x = reorder(Group.1, -mean), y = mean, fill = Algae_group)) +     
  geom_bar(stat = "identity", colour="black")+
  scale_fill_manual(values=as.character(colors$color)) +
  labs(y=paste("Average relative abundance(%)"), x=paste("Phylum"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())
p

##supplementary figure 4

level <- which(colnames(taxonomy)=="Genus_clean")
asv_table_fam<- aggregate(asv_table_norm, by=list(taxonomy[, level]), FUN=sum)
otu_long <- melt(asv_table_fam, id.vars = "Group.1", variable.name = "Sample")
points <- cbind(otu_long, design[match(otu_long$Sample, design$SampleID), ])
points2 <- points[points$Group.1 %in% c("Chlamydomonas"), ]

sd<- ddply(points2, c("Group.1", "Site"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

p<- ggplot(sd, aes(x = reorder(Site, -mean), y = mean)) +     
  geom_bar(stat = "identity", colour="black")+
  #scale_fill_manual(values=as.character(colors$color)) +
  labs(y=paste("Average relative abundance of Chlamydomonas sp.(%)"), x=paste("SW sites"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.9))+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())
p
