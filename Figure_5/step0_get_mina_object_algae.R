
# Load data and get 'mina' object
## 2025.09.03
#cls <- read.table("20250818_clusters_env_factors.txt", header=T, sep="\t")
## in latested cluster file, only 100 sites are included; while previous design
## file has 150 sites

asv <- read.table("asv_core_plants.txt", header=T, sep="\t")
design <- read.table("merged_design.txt", header=T, sep="\t")

cls2 <- unique(cls[, c("Site.1", "Cluster")])
## filter out NA in column "Cluster" for bs_pm()
design <- design[!is.na(design$Cluster), ]

# Filter ASVs present in more than 5% of samples
# keep_asv <- rowSums(asv > 0) > (0.05 * ncol(asv))
# asv <- asv[keep_asv, ]
# 615 410

bac_tax <- read.delim2("taxonomy_mod_bac2.txt", header=T, sep="\t")
fun_tax <- read.table("taxonomy_mod_fun2.txt", header=T, sep="\t")

bac_tax <- bac_tax[bac_tax$FeatureID %in% rownames(asv), ]
fun_tax <- fun_tax[fun_tax$Feature.ID %in% rownames(asv), ]

# write.table(bac_tax, "taxonomy_mod_bac2_filtered.txt", 
#             sep="\t", row.names=F, quote=F)
# write.table(fun_tax, "taxonomy_mod_fun2_filtered.txt", 
#             sep="\t", row.names=F, quote=F)

## After filtering, all left fungal ASVs are Core.
colnames(fun_tax)[1] <- "FeatureID"
colnames(bac_tax)[2] <- "Kingdom"
bac_tax_simp <- bac_tax[, c("FeatureID", "Kingdom", "Order")]
fun_tax_simp <- fun_tax[, c("FeatureID","Kingdom","Order")]

tax <- rbind(bac_tax_simp, fun_tax_simp)

# Create new mina object
algae <- new("mina", tab = as.matrix(asv), des = design)
algae <- norm_tab(algae, method = "raref")
algae <- fit_tabs(algae)

algae <- get_rep(algae)