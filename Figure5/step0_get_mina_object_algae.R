# Load data and get 'mina' object
## 2025.09.03
cls <- read.table("../data/20250818_clusters_env_factors.txt", header=T, sep="\t")
## in latest cluster file, only 100 sites are included; while previous design
## file has 150 sites

asv <- read.table("../data/asv_core_plants.txt", header=T, sep="\t")
design <- read.table("../data/merged_design.txt", header=T, sep="\t")

cls2 <- unique(cls[, c("Site.1", "Cluster")])

design <- design[design$Site %in% cls2$Site.1, ]
## 303 samples in 100 sites

asv1 <- asv[, colnames(asv) %in% design$Sample_ID]
dim(asv)
## 25515 410
dim(asv1)

##################
## 25515 280 @Paloma: is this correct? 280 samples out of 303 samples in ASV table
## examples of samples filtered out here:
## "SurveyECOGEN2023_201" "SurveyECOGEN2023_249" "SurveyECOGEN2023_297" 
## "SurveyECOGEN2023_345" "SurveyECOGEN2023_209"


## Ignore the "Cluster" column in design, use the one in 
## 20250818_clusters_env_factors.txt

design <- design[, colnames(design) != "Cluster"]
design <- merge(design, cls2, by.x="Site", by.y="Site.1", all.x=T)
dim(design)
## 303 6

# Filter ASVs present in more than 5% of samples
# keep_asv <- rowSums(asv > 0) > (0.05 * ncol(asv))
# asv <- asv[keep_asv, ]
# 615 410

bac_tax <- read.delim2("../data/taxonomy_mod_bac2.txt", header=T, sep="\t")
fun_tax <- read.delim2("../data/taxonomy_mod_fun2.txt", header=T, sep="\t")

bac_tax <- bac_tax[bac_tax$FeatureID %in% rownames(asv1), ]
fun_tax <- fun_tax[fun_tax$Feature.ID %in% rownames(asv1), ]

## only keep bacterial and fungal ASVs that are in "core"
core_bac <- c("Burkholderiales", "Caulobacterales", "Microtrichales",
              "Propionibacteriales", "Rhizobiales", "Solirubrobacterales",
              "Sphingomonadales")

bac_tax_core <- bac_tax[bac_tax$Order %in% core_bac, ]

core_fun <- c("Chaetothyriales", "Helotiales", "Pleosporales",
              "Hypocreales", "Xylariales")
fun_tax_core <- fun_tax[fun_tax$Order %in% core_fun, ]

asv_core <- asv1[rownames(asv1) %in% bac_tax_core$FeatureID | 
             rownames(asv1) %in% fun_tax_core$Feature.ID, ]
dim(asv_core)
## 19649   280

## After filtering, all left fungal ASVs are Core.
colnames(fun_tax_core)[1] <- "FeatureID"
colnames(bac_tax_core)[2] <- "Kingdom"
bac_tax_simp <- bac_tax_core[, c("FeatureID", "Kingdom", "Order", "Family")]
fun_tax_simp <- fun_tax_core[, c("FeatureID","Kingdom","Order", "Family")]

tax <- rbind(bac_tax_simp, fun_tax_simp)

### ASVs with Family tax assignments
tax_fam <- tax[grepl("f__", tax$Family, ignore.case = T) == T, ]
dim(tax_fam)
## 18377   4

asv_core <- asv_core[rownames(asv_core) %in% tax_fam$FeatureID, ]
dim(asv_core)
## 18377   280


# Filter core ASVs present in less than 5% of samples
keep_asv <- rowSums(asv_core > 0) > (0.05 * ncol(asv_core))
asv_core_keep <- asv_core[keep_asv, ]
dim(asv_core_keep)

# 490 280
library(ggplot2)
t <- data.frame(Sample_ID = colnames(asv_core_keep), 
                Reads = colSums(asv_core_keep))

## reads count distribution
ggplot(t, aes(x=Reads)) + geom_histogram(binwidth = 500) +
  theme_bw() + geom_vline(xintercept = 1000, color="salmon", linetype = "dashed") 

################################################################################
# Create new mina object
library(mina)
algae <- new("mina", tab = as.matrix(asv_core_keep), des = design)
algae <- norm_tab(algae, method = "raref", depth = 500)

algae <- fit_tabs(algae)

#algae <- get_rep(algae)
