## Supplementary figure 2 
# Every plots can be generated individually 
# Packages

library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggdendro)
library(cowplot)
library(GUniFrac)
library(paletteer)
library(scales)
library(patchwork)

# CPCoA functions

library(vegan)
#~ library(calibrate)
#~ library(Biostrings)
#~ library(ape)

variability_table <- function(cca){
  
  chi <- c(cca$tot.chi,
           cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
  
}

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig/eig_tot
  return(var_propdf)
}

pca_var_props <- function(cca){
  
  eig_tot <- sum(cca$CA$eig)
  var_propdf <- cca$CA$eig/eig_tot
  return(var_propdf)
}

cca_ci <- function(cca, permutations=5000){
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
  return(ci)
  
}

#sources(cpcoa.func.R)


# Supplementary figure Xa) CPCoA of algal community 

asv_18S=read.delim("ASV_18S.tsv", row.names=1, check.names=FALSE) #load data
taxo_18S=read.delim("TAXO_18S.tsv")
map_18S=read.delim("mapping.txt")

mat_18S=asv_18S[, colnames(asv_18S) %in% map_18S$SampleID] #filtering and normalization
mat_18S=mat_18S[, colSums(mat_18S) > 0]
asv_table_norm_18S=apply(mat_18S, 2, function(x) x/sum(x))

idx=map_18S$SampleID %in% colnames(asv_table_norm_18S) #link tables
map_18S=map_18S[idx, ]
map_18S=map_18S[match(colnames(asv_table_norm_18S), map_18S$SampleID), ]
bray_curtis_18S=vegdist(t(asv_table_norm_18S), method="bray")#bray-curtis distance for the PCoA

capscale_18S <- capscale(bray_curtis_18S ~ Location+Soil, data=map_18S, add=F, sqrt.dist=T) #constrained PCoA and statistics
perm_anova_18S=anova.cca(capscale_18S)
var_tbl18S <- variability_table(capscale_18S)
variance <- var_tbl18S["constrained", "proportion"]
p.val <- perm_anova_18S[1, 4]
eig=capscale_18S$CCA$eig
p_val=perm_anova_18S[1, 4]

points <- capscale_18S$CCA$wa[, 1:2] #extract coordinates
points <- as.data.frame(points)
colnames(points)=c("x", "y")
points$SampleID=rownames(points)
points=left_join(points, map_18S, by="SampleID")#join metadata

centroids=aggregate(cbind(x, y) ~ Location+Soil, data=points, FUN=mean)#calculate centroids and segments
names(centroids)=c("Location","Soil", "cen_x", "cen_y")
segments=left_join(points, centroids, by=c("Location","Soil"))

#plot
cpcoa_18S=ggplot(segments, aes(x=x, y=y, color=Location,shape=Soil)) +
  geom_segment(aes(xend=cen_x, yend=cen_y), alpha=0.5) +
  geom_point(size=3) +
  geom_point(data=centroids, aes(x=cen_x, y=cen_y), shape=3, size=5, stroke=1.5, color="black") +
  scale_color_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) + 
  labs(x=paste0("cPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)"),
       y=paste0("cPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)"))+
  ggtitle(paste(format(100 * variance, digits=2), " % of variance; p=",format(p.val, digits=2),sep="")) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), legend.text=element_text(size=12), 
        legend.title=element_text(size=14), axis.title=element_text(size=14),)


# Supplementary figure 2b) Algae relative abundance

asv=read.table("ASV_18S.tsv", row.names=1, check.names=FALSE, header=TRUE)
taxo=read.delim("TAXO_18S.tsv")
map=read.delim("mapping.txt")

colnames(taxo)[1]="FeatureID"
taxo=separate(taxo, Taxon, c("D","Phylum","C","O","F","G","S"), sep=";", fill="right")#correct dataframe structure
taxo$Phylum=trimws(gsub("p__", "", taxo$Phylum ))#remove p__ before phylum names
asv_rel=apply(asv, 2, function(x) x/sum(x))#relative abundance 

algae=c("Chlorophyta", "Charophyta")#select algae only
est_algue=taxo$Phylum[match(rownames(asv_rel), taxo$FeatureID)] %in% algae
algue_sum=colSums(asv_rel[est_algue, , drop=FALSE], na.rm=TRUE)

tab=data.frame(SampleID=names(algue_sum), Abundance=algue_sum)
tab=left_join(tab, map[, c("SampleID", "Location", "Soil")], by="SampleID")#keep Location and Soil
tab=filter(tab, !is.na(Location), !is.na(Soil))#remove NA

#plot
violinplot=ggplot(tab, aes(x=Location, y=Abundance, fill=Location)) +
  geom_violin(alpha=0.4) +  
  geom_point(aes(shape=Soil, fill=Location), color="black", size=2.5, position=position_jitter(width=0.1)) +  
  annotate("text", x=1.5, y=0.12, label="ns", size=4) +
  scale_fill_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) +
  scale_shape_manual(values=c(21, 24, 22)) +
  labs(y="Relative abundance (%)", x="") +
  theme(axis.text.x=element_text(size=12, colour="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14), 
        axis.title=element_text(size=14))


# Supplementary figure 2c) CPCoA of bacterial community

#same as a) and e) 
asv_16S=read.delim("ASV_16S.tsv", row.names=1, check.names=FALSE)
taxo_16S=read.delim("TAXO_16S.tsv")
map_16S=read.delim("mapping.txt")
mat_16S=asv_16S[, colnames(asv_16S) %in% map_16S$SampleID]
mat_16S=mat_16S[, colSums(mat_16S) > 0]
asv_table_norm_16S=apply(mat_16S, 2, function(x) x/sum(x))

idx=map_16S$SampleID %in% colnames(asv_table_norm_16S)
map_16S=map_16S[idx, ]
map_16S=map_16S[match(colnames(asv_table_norm_16S), map_16S$SampleID), ]
bray_curtis_16S=vegdist(t(asv_table_norm_16S), method="bray")

capscale_16S <- capscale(bray_curtis_16S ~ Location+Soil, data=map_16S, add=F, sqrt.dist=T)
perm_anova_16S=anova.cca(capscale_16S)
var_tbl16S <- variability_table(capscale_16S)
variance <- var_tbl16S["constrained", "proportion"]
p.val <- perm_anova_16S[1, 4]
eig=capscale_16S$CCA$eig
p_val=perm_anova_16S[1, 4]

points <- capscale_16S$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points)=c("x", "y")


points$SampleID=rownames(points)
points=left_join(points, map_16S, by="SampleID")
centroids=aggregate(cbind(x, y) ~ Location+Soil, data=points, FUN=mean)
names(centroids)=c("Location","Soil", "cen_x", "cen_y")
segments=left_join(points, centroids, by=c("Location","Soil"))


cpcoa16S=ggplot(segments, aes(x=x, y=y, color=Location,shape=Soil)) +
  geom_segment(aes(xend=cen_x, yend=cen_y), alpha=0.5) +
  geom_point(size=3) +
  geom_point(data=centroids, aes(x=cen_x, y=cen_y), shape=3, size=5, stroke=1.5, color="black") +
  scale_color_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) + 
  labs(x=paste0("cPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)"),
       y=paste0("cPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)"))+
  ggtitle(paste(format(100 * variance, digits=2), " % of variance; p=",format(p.val, digits=2),sep="")) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), legend.text=element_text(size=12), 
        legend.title=element_text(size=14), axis.title=element_text(size=14),)


# Supplementary figure 2d) Relative abundance of bacterial core orders

asv_16S=read.delim("ASV_16S.tsv", row.names=1, check.names=FALSE)
taxo_16S=read.delim("TAXO_16S.tsv")
map16S=read.delim("mapping.txt")

asv_16S=asv_16S[, colSums(asv_16S) >= 1000]#remove everything under 1000 

colnames(taxo_16S)[1]="FeatureID"
taxo_16S=separate(taxo_16S, Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", fill="right")
taxo_16S$Order=trimws(gsub("o__", "", taxo_16S$Order))
taxo_16S$Order[taxo_16S$Order=="Hyphomicrobiales"] <- "Rhizobiales" #fix Rhizobiales name 

coregrp16S=c("Burkholderiales", "Caulobacterales", "Microtrichales", "Propionibacteriales", "Rhizobiales", "Solirubrobacterales", "Sphingomonadales")#select known core members
taxo_core16S=taxo_16S[taxo_16S$Order %in% coregrp16S, ]
mat_core16S=asv_16S[rownames(asv_16S) %in% taxo_core16S$FeatureID, ]

asv_rel16S=apply(mat_core16S, 2, function(x) x/sum(x))#relative abundance
asv_long16S=as.data.frame(asv_rel16S)
asv_long16S$FeatureID=rownames(asv_long16S)
asv_long16S=pivot_longer(asv_long16S, -FeatureID, names_to="SampleID", values_to="Abundance")

barplot_data16S=left_join(asv_long16S, taxo_core16S[, c("FeatureID", "Order")], by="FeatureID")#barplot dataset
barplot_data16S=left_join(barplot_data16S, map16S[, c("SampleID", "Location")], by="SampleID")
barplot_data16S=barplot_data16S %>% filter(!is.na(Order)) %>% filter(!is.na(Location)) %>% filter(Order != "NA")

data_mean=aggregate(Abundance ~ Location + Order, data=barplot_data16S, FUN=mean)#mean calculation
data_sd=aggregate(Abundance ~ Location + Order, data=barplot_data16S, FUN=sd)#standard deviation
data_n=aggregate(Abundance ~ Location + Order, data=barplot_data16S, FUN=length)#data lenght

data_final=data.frame(data_mean, SD=data_sd$Abundance, N=data_n$Abundance)
data_final$SE=data_final$SD / sqrt(data_final$N)#final dataset
annot_df=aggregate(Abundance + SE ~ Order, data=data_final, max)
colnames(annot_df)[2]="y_pos"

annot_df$Label=c("***", "***", "***", "***", "***", "***", "***")#wilcoxon results annotation

#plot
barplot16S=ggplot(data_final, aes(x=Order, y=Abundance, fill=Location)) +
  geom_bar(stat="identity", position="dodge", col="black") +
  geom_errorbar(aes(ymin=Abundance, ymax=Abundance+SE), width=0.2, position=position_dodge(width=0.8)) +
  geom_text(data=annot_df, aes(x=Order, y=y_pos, label=Label), inherit.aes=FALSE, size=6) +
  scale_fill_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) +
  theme(axis.text.x=element_text(size=12, colour="black", angle=45, hjust=1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14), 
        axis.title=element_text(size=14))+ylab("Average relative abundance (%)")

# Supplementary figure 2f) CPCoA of fungal community

#same as a) and c)

asv_ITS=read.delim("ASV_ITS.tsv", row.names=1, check.names=FALSE)
taxo_ITS=read.delim("TAXO_ITS.tsv")
map_ITS=read.delim("mapping.txt")
mat_ITS=asv_ITS[, colnames(asv_ITS) %in% map_ITS$SampleID]
mat_ITS=mat_ITS[, colSums(mat_ITS) > 0]
asv_table_norm_ITS=apply(mat_ITS, 2, function(x) x/sum(x))

idx=map_ITS$SampleID %in% colnames(asv_table_norm_ITS)
map_ITS=map_ITS[idx, ]
map_ITS=map_ITS[match(colnames(asv_table_norm_ITS), map_ITS$SampleID), ]
bray_curtis_ITS=vegdist(t(asv_table_norm_ITS), method="bray")

capscale_ITS <- capscale(bray_curtis_ITS ~ Location+Soil, data=map_ITS, add=F, sqrt.dist=T)
perm_anova_ITS=anova.cca(capscale_ITS)
print(perm_anova_ITS)
var_tblITS <- variability_table(capscale_ITS)
variance <- var_tblITS["constrained", "proportion"]
p.val <- perm_anova_ITS[1, 4]
eig=capscale_ITS$CCA$eig
p_val=perm_anova_ITS[1, 4]

points <- capscale_ITS$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points)=c("x", "y")
points$SampleID=rownames(points)
points=left_join(points, map_ITS, by="SampleID")

centroids=aggregate(cbind(x, y) ~ Location+Soil, data=points, FUN=mean)
names(centroids)=c("Location","Soil", "cen_x", "cen_y")
segments=left_join(points, centroids, by=c("Location","Soil"))

cpcoaITS=ggplot(segments, aes(x=x, y=y, color=Location,shape=Soil)) +
  geom_segment(aes(xend=cen_x, yend=cen_y), alpha=0.5) +
  geom_point(size=3) +
  geom_point(data=centroids, aes(x=cen_x, y=cen_y), shape=3, size=5, stroke=1.5, color="black") +
  scale_color_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) + 
  labs(x=paste0("cPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)"),
       y=paste0("cPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)"))+
  ggtitle(paste(format(100 * variance, digits=2), " % of variance; p=",format(p.val, digits=2),sep="")) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), legend.text=element_text(size=12), 
        legend.title=element_text(size=14), axis.title=element_text(size=14),)

# Supplementary figure 2g) Relative abundance of fungal core members

#same as d)
asv_ITS=read.delim("ASV_ITS.tsv", row.names=1, check.names=FALSE)
taxo_ITS=read.delim("TAXO_ITS.tsv")
mapITS=read.delim("mapping.txt")
asv_ITS=asv_ITS[, colSums(asv_ITS)>=1000]

colnames(taxo_ITS)[1]="FeatureID"
taxo_ITS=separate(taxo_ITS, Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", fill="right")
taxo_ITS$Order=trimws(gsub("o__", "", taxo_ITS$Order))

coregrpITS=c("Chaetothyriales", "Helotiales", "Hypocreales", "Pleosporales", "Xylariales")
taxo_coreITS=taxo_ITS[taxo_ITS$Order %in% coregrpITS, ]
mat_coreITS=asv_ITS[rownames(asv_ITS) %in% taxo_coreITS$FeatureID, ]

asv_relITS=apply(mat_coreITS, 2, function(x) x/sum(x))
asv_longITS=as.data.frame(asv_relITS)
asv_longITS$FeatureID=rownames(asv_longITS)
asv_longITS=pivot_longer(asv_longITS, -FeatureID, names_to="SampleID", values_to="Abundance")

barplot_dataITS=left_join(asv_longITS, taxo_coreITS[, c("FeatureID", "Order")], by="FeatureID")
barplot_dataITS=left_join(barplot_dataITS, mapITS[, c("SampleID", "Location")], by="SampleID")
barplot_dataITS=barplot_dataITS %>% filter(!is.na(Order)) %>% filter(!is.na(Location)) %>% filter(Order != "NA")

data_mean=aggregate(Abundance ~ Location + Order, data=barplot_dataITS, FUN=mean)
data_sd=aggregate(Abundance ~ Location + Order, data=barplot_dataITS, FUN=sd)
data_n=aggregate(Abundance ~ Location + Order, data=barplot_dataITS, FUN=length)

data_final=data.frame(data_mean, SD=data_sd$Abundance, N=data_n$Abundance)
data_final$SE=data_final$SD / sqrt(data_final$N)

annot_df=aggregate(Abundance + SE ~ Order, data=data_final, max)
colnames(annot_df)[2]="y_pos"

annot_df$Label=c("ns", "ns", "ns", "ns", "*")

barplot_ITS=ggplot(data_final, aes(x=Order, y=Abundance, fill=Location)) +
  geom_bar(stat="identity", position="dodge", col="black") +
  geom_errorbar(aes(ymin=Abundance-SE, ymax=Abundance+SE), width=0.2, position=position_dodge(width=0.8)) +
  geom_text(data=annot_df, aes(x=Order, y=y_pos, label=Label), inherit.aes=FALSE, size=6) +
  scale_fill_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) +
  theme(axis.text.x=element_text(size=12, colour="black", angle=45, hjust=1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14), 
        axis.title=element_text(size=14))+ylab("Average relative abundance (%)")

# Supplementary figure 2h) Aggregated relative abundance of bacterial and fungal core and non-core communities 

#analysis for bacteria
asv16S=read.delim("ASV_16S.tsv", row.names=1, check.names=FALSE)
taxo16S=read.delim("TAXO_16S.tsv")
map16S=read.delim("mapping.txt")

taxo16S=separate(taxo16S, Taxon, c("D","P","C","Order"), sep=";", fill="right")
taxo16S$Order=gsub("o__", "", taxo16S$Order)
taxo16S$Order=trimws(taxo16S$Order)

coregrp16S=c("Burkholderiales", "Caulobacterales", "Microtrichales", "Propionibacteriales", "Hyphomicrobiales", "Solirubrobacterales", "Sphingomonadales")#select core members

asv_rel16S=apply(asv16S, 2, function(x) x/sum(x))
asv_df16S=as.data.frame(asv_rel16S)
asv_df16S$FeatureID=rownames(asv_df16S)

df_long16S=pivot_longer(asv_df16S, -FeatureID, names_to="SampleID", values_to="Abundance")

df16S=left_join(df_long16S, taxo16S[, c("FeatureID", "Order")], by="FeatureID")
df16S=left_join(df16S, map16S[, c("SampleID", "Location")], by="SampleID")

df16S$Group="Non-Core"
df16S$Group[df16S$Order %in% coregrp16S]="Core"

final16S=aggregate(Abundance ~ SampleID + Location + Group, data=df16S, sum)

#analysis for fungi
asvITS=read.delim("ASV_ITS.tsv", row.names=1, check.names=FALSE)
taxoITS=read.delim("TAXO_ITS.tsv")
mapITS=read.delim("mapping.txt")

taxoITS=separate(taxoITS, Taxon, c("D","P","C","Order"), sep=";", fill="right")
taxoITS$Order=gsub("o__", "", taxoITS$Order)
taxoITS$Order=trimws(taxoITS$Order)

coregrpITS=c("Chaetothyriales", "Helotiales", "Hypocreales", "Pleosporales", "Xylariales")#select core members

asv_relITS=apply(asvITS, 2, function(x) x/sum(x))
asv_dfITS=as.data.frame(asv_relITS)
asv_dfITS$FeatureID=rownames(asv_dfITS)

df_longITS=pivot_longer(asv_dfITS, -FeatureID, names_to="SampleID", values_to="Abundance")

dfITS=left_join(df_longITS, taxoITS[, c("FeatureID", "Order")], by="FeatureID")
dfITS=left_join(dfITS, mapITS[, c("SampleID", "Location")], by="SampleID")

dfITS$Group="Non-Core"
dfITS$Group[dfITS$Order %in% coregrpITS]="Core"

finalITS=aggregate(Abundance ~ SampleID + Location + Group, data=dfITS, sum)
final16S$Label=c("***", "***")

#plots
p16S=ggplot(final16S, aes(x=Group, y=Abundance, fill=Location)) +
  geom_boxplot(show.legend=FALSE) +
  annotate("text", x=1, y=0.4 , label="*", size=8) +
  annotate("text", x=2, y=1 , label="*", size=8) +
  scale_fill_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) +
  theme(axis.text.x=element_text(size=12, colour="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14), 
        axis.title=element_text(size=14))+
  ggtitle("Bacteria")+ylab("Relative abundance (%)")


pITS=ggplot(finalITS, aes(x=Group, y=Abundance, fill=Location)) +
  geom_boxplot() +
  annotate("text", x=1, y=0.7 , label="ns", size=5) +
  annotate("text", x=2, y=0.8 , label="ns", size=5) +
  scale_fill_manual(values=c("Deep"="#666889", "Surface"="#6caa91")) +
  theme(axis.text.x=element_text(size=12, colour="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14), 
        axis.title=element_text(size=14))+
  ggtitle("Fungi")+ylab("Relative abundance (%)")


boxplot=plot_grid(p16S, pITS)


# Statistical analysis (Wilcoxon Test , Holm-corrected)

# Core vs Non-Core (bacteria and fungi)

# Bacteria 16S
final16S$LocGroup <- paste(final16S$Location, final16S$Group, sep="_")
pairwise.wilcox.test(final16S$Abundance, final16S$LocGroup)

# Fungi ITS
finalITS$LocGroup <- paste(finalITS$Location, finalITS$Group, sep="_")
pairwise.wilcox.test(finalITS$Abundance, finalITS$LocGroup)


# Individual orders by location (fungi)
print("Individual orders by location (fungi)")
df_its <- aggregate(Abundance ~ SampleID + Location + Order, data=barplot_dataITS, sum)

print(wilcox.test(Abundance ~ Location, data=subset(df_its, Order == "Chaetothyriales")))
print(wilcox.test(Abundance ~ Location, data=subset(df_its, Order == "Helotiales")))
print(wilcox.test(Abundance ~ Location, data=subset(df_its, Order == "Hypocreales")))
print(wilcox.test(Abundance ~ Location, data=subset(df_its, Order == "Pleosporales")))
print(wilcox.test(Abundance ~ Location, data=subset(df_its, Order == "Xylariales")))


# Individual orders by location (bacteria)
print("Individual orders by location (bacteria)")
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Burkholderiales")))
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Caulobacterales")))
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Hyphomicrobiales")))
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Microtrichales")))
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Propionibacteriales")))
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Solirubrobacterales")))
print(wilcox.test(Abundance ~ Location, data=subset(df16S, Order == "Sphingomonadales")))



#all plots

final_plot <- (
  (cpcoa_18S | violinplot) / 
    (cpcoa16S  | barplot16S) / 
    (cpcoaITS  | barplot_ITS) / 
    (plot_spacer() + boxplot + plot_spacer() + plot_layout(widths = c(0.5, 2, 0.5)))) + plot_annotation(tag_levels = 'a') 

final_plot

ggsave(final_plot,filename="supplementary_figure_2.pdf",path="H:/Priv/PROJET_MICROBIOTA/Scripts/SuppFigX" ,width = 20, height = 20)
