##supplementary figure 5
##data in figure 1 folder

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

design<- read.table("design_noempty.txt", header=TRUE, sep="\t")
asv_table<- read.table("otu_table_algae.txt", header=TRUE, sep="\t")

asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

design_cluster <- design[design$environmental_cluster %in% c("1", "2", "3"), ]
idx <- colnames(asv_table_norm) %in% design_cluster$SampleID

bray_curtis <- vegdist(t(asv_table_norm[,idx]), method="bray")
bray_curtis<- as.matrix(bray_curtis)

idx <-design_cluster$SampleID %in% colnames(bray_curtis)
design_cluster<- design_cluster[idx,]


source("cpcoa.func.R")
sqrt_transform <- T
capscale.rate <- capscale(bray_curtis ~ environmental_cluster,
                          data=design_cluster, add=F, sqrt.dist=sqrt_transform)
# ANOVA-like permutation analysis

perm_anova.rate <- anova.cca(capscale.rate)
print(perm_anova.rate)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.rate <- variability_table(capscale.rate)

eig <- capscale.rate$CCA$eig

variance <- var_tbl.rate["constrained", "proportion"]
p.val <- perm_anova.rate[1, 4]

# extract the weighted average (sample) scores

points <- capscale.rate$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_cluster[match(rownames(points), design_cluster$SampleID), ])

colors <- data.frame(group=c("1","2", "3"),
                     color=c("#fb8072","#6AA84F", "#47A5EC"))

p <- ggplot(points, aes(x=x, y=y, color=environmental_cluster)) +
  geom_point(alpha=.7, size=6)+
  scale_color_manual(values=colors$color)+
  labs(x=paste("cPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("cPCoA 2 (", format(100 * eig[2] / sum(eig),	 digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=2), " % of variance; p=",
                format(p.val, digits=2),
                sep="")) +
  theme_classic()
p
