# Network Construction

library(pheatmap)   
library(Hmisc)      
library(stringr)    
library(dplyr)      
library(tibble)     
library(reshape2)
library(vegan)


taxonomy_bac <- read.table("taxonomy_mod_bac.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)
taxonomy_fun <- read.table("taxonomy_mod_fun.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)

taxonomy_bac  <- taxonomy_bac  %>% rename(Kingdom = Domain)
cols_to_keep <- c("Feature.ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy_fun <- taxonomy_fun[, cols_to_keep]
taxonomy_bac <- taxonomy_bac[, cols_to_keep]

taxonomy_all <- rbind(taxonomy_bac , taxonomy_fun)
# Looks like that : Kingdom | Phylum | Class | Order | Family | Genus
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")  

tax_cols <- tax_cols[tax_cols %in% colnames(taxonomy_all)]  

taxonomy_all <- taxonomy_all %>%
  mutate(Taxonomy = paste(!!!syms(tax_cols), sep=" | ")) %>%
  select(Feature.ID, Taxonomy)  

# Modify for each bootstrap:
# asv_norm_raref --> boot1-20
merged_data <- left_join(ASVs_cluster3, taxonomy_all, by = "Feature.ID")

merged_data <- na.omit(merged_data)

merged_data$Taxonomy <- make.unique(as.character(merged_data$Taxonomy))

rownames(merged_data) <- merged_data$Taxonomy

merged_data <- merged_data %>% select(-Taxonomy)

merged_data <- merged_data %>% select(-Feature.ID)

merged_data <- as.matrix(merged_data)

taxa_merged <- str_split_fixed(rownames(merged_data), "\\|", 6)
colnames(taxa_merged) <- c("Kindom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(taxa_merged) <- rownames(merged_data)



merged_data <-  as.matrix(apply(merged_data, 2, function(x) x/sum(x)))

n_sample <- ncol(merged_data)
n_sample_filter <- round(n_sample/20)

ra_mat_filter1 <- merged_data[rowSums(merged_data > 0) > n_sample_filter, ]
ra_mat_filter2 <- merged_data[rowSums(merged_data) / n_sample > 0.00001, ]

pass_filter <- intersect(rownames(ra_mat_filter1), rownames(ra_mat_filter2))
rm(ra_mat_filter1, ra_mat_filter2)

ra_mat_filter <- merged_data[pass_filter, ]

# Re-normalize the relative abundance by column
ra_mat_filter <-  as.matrix(apply(ra_mat_filter, 2, function(x) x/sum(x)))

ra_corr <- rcorr(t(ra_mat_filter), type = "spearman")

# Correlation matrix
mat_cor <- ra_corr$r

# Significance matrix
sig_cor <- ra_corr$P

# Filter out the medium/strong (abs > 0.2) and significant (P < 0.01) correlation
mat_cor[abs(mat_cor) < 0.2] <- 0
mat_cor[sig_cor > 0.01] <- 0
diag(mat_cor) <- 0

mat_cor_filter <- mat_cor[rowSums(mat_cor) > 0, colSums(mat_cor) > 0]
sig_cor_filter <- sig_cor[rowSums(mat_cor) > 0, colSums(mat_cor) > 0]

mat_cor_filter2 <- mat_cor_filter



this_taxa <- taxa_merged[rownames(taxa_merged) %in% rownames(mat_cor_filter), ]
this_taxa_color <- this_taxa[, 2]

taxa_map <- data.frame(Full_name = rownames(this_taxa),
                       Kingdom = this_taxa[, 2],
                       Order = this_taxa[, 4])
taxa_map <- taxa_map[match(taxa_map$Full_name, rownames(mat_cor_filter)), ]
rownames(mat_cor_filter) <- colnames(mat_cor_filter) <- taxa_map$Order

taxa_map2 <- data.frame(Kingdom = factor(x = taxa_map$Kingdom,
                                         levels = unique(taxa_map$Kingdom)))
rownames(taxa_map2) <- make.unique(as.character(taxa_map$Order))

taxa_map2$index <- rownames(taxa_map2)

taxa_map2 <- taxa_map2[!grepl("\\.\\d+$", taxa_map2$index), ]

rownames(taxa_map2) <- taxa_map2$index
taxa_map2$index <- NULL
mat_cor_filter2[upper.tri(mat_cor_filter2)] <- 0

# Create the node and the edge
df1 <- melt(mat_cor_filter2, na.rm = TRUE)
colnames(df1) <- c("node1", "node2", "edge")

df2 <- melt(sig_cor_filter, na.rm = TRUE)
colnames(df2) <- c("node1", "node2", "sig")

df_net <- merge(df1, df2, by = c("node1", "node2"))
df_net_sig <- df_net[df_net$sig < 0.01, ]
df_net_sig <- df_net_sig[df_net_sig$edge != 0, ]

this_taxa <- taxa_merged[rownames(taxa_merged) %in% rownames(mat_cor_filter2), ]
this_taxa_color <- this_taxa[, 2]

taxa_map <- data.frame(Full_name = rownames(this_taxa),
                       Kingdom = this_taxa[, 1],
                       Phylum = this_taxa[, 2],
                       Class = this_taxa[, 3],
                       Order = this_taxa[, 4], 
                       Family = this_taxa[,5],
                       Genus = this_taxa[,6])
taxa_map <- taxa_map[match(taxa_map$Full_name, rownames(mat_cor_filter2)), ]
dim(taxa_map)
dim(df_net_sig)



# Download
write.table(taxa_map, "node_C3.txt", quote = F, sep = "\t", row.names = F)
write.table(df_net_sig, "edge_C3.txt", quote = F, sep = "\t", row.names = F)
