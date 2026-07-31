library(mina)
library(RSpectra)
library(reshape2)
library(stringr)
library(vegan)
library(ggplot2)

source("20250801_step0_get_mina_object_algae.R")
## 570 336 

if (!dir.exists("./bs_pm_out2/")) dir.create("./bs_pm_out2/")
    
algae <- bs_pm(algae, group = "Cluster", rm = F, individual = T, 
               out_dir = "./bs_pm_out2/",
               g_size = 30, s_size = 15, bs = 11, pm = 11)

################################################################################
### For feature selection
if (!dir.exists("./bs_pm_out2_dis/")) dir.create("./bs_pm_out2_dis/")
algae2 <- net_dis_indi("./bs_pm_out2/", method = "spectra",
                       dir = "./bs_pm_out2_dis/")

dis_stat(algae2)

# Compare Distance_Mean Distance_SD Distance_PM_Mean Distance_PM_SD     N          p
# 1     1_1      12.38969    4.242427         12.53519       4.187272 53361 0.51504816
# 2     1_2      13.88314    4.444489         13.47354       3.840211 14641 0.47889633
# 3     1_3      23.03685    6.156791         12.03207       3.471797 14641 0.05825707
# 4     2_2      14.01978    4.616577         11.80612       3.300598 53361 0.36124958
# 5     2_3      22.75499    6.116734         12.71419       3.289231 14641 0.06214998
# 6     3_3      14.57846    4.533543         13.01852       4.141338 53361 0.40375173


### Panel a, PCA of bootstrap network
################################################################################

bs_files <- list.files("./bs_pm_out2_dis/",
                       pattern = "^spectra_bs_", full.names = TRUE)

## for bootstrap results and plot
bs_vectors <- c()
for (bs in bs_files) {
  this <- readRDS(bs)
  bs_vectors <- rbind(bs_vectors, this)
}

dis_bs <- as.matrix(vegdist(bs_vectors, method = "euclidean"))
rownames(dis_bs) <- colnames(dis_bs) <- paste0(rownames(dis_bs), "_", seq(1 : nrow(dis_bs)))

## get the group information
idx <- seq(1,length(colnames(dis_bs)) * 2, by = 2)
design <- data.frame(Net = colnames(dis_bs),
                     Group = unlist(strsplit(colnames(dis_bs), "_bs"))[idx])

ad <- adonis2(formula = dis_bs ~ Group,
              data = design, by = "margin", add = F, parallel = 40)
print(ad)

# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dis_bs ~ Group, data = design, add = F, by = "margin", parallel = 40)
# Df SumOfSqs      R2      F Pr(>F)    
# Group      2    14948 0.54982 78.776  0.001 ***
#   Residual 129    12239 0.45018                  
# Total    131    27188 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


dis_bs_dmr <- cmdscale(dis_bs, k = 4, eig = T)
eig <- dis_bs_dmr$eig
eig[eig < 0] <- 0
eig1 <- 100 * eig[1] / sum(eig)
eig2 <- 100 * eig[2] / sum(eig)

points <- as.data.frame(dis_bs_dmr$points[, 1:2])
colnames(points) <- c("x", "y")
points$Net <- rownames(points)
points <- merge(points, design)

p_a <- ggplot(points, aes(x, y, color = factor(Group))) + #, shape = factor(Group))) +
  geom_point(size = 3, alpha = 0.8, shape = 16) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("salmon",  "olivedrab", "deepskyblue")) +
  # scale_shape_manual(values = c(15, 16, 17, 18),  labels = paste("Cluster", 1:4)) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top",legend.title = element_blank()) +
  labs(x = paste0("PC1 (", format(eig1, digits = 4), "%)"), 
       y = paste0("PC2 (", format(eig2, digits = 4), "%)"))

p_a
ggsave("Figure_5a.pdf", width = 5, height = 3)
################################################################################


################################################################################
## modified from get_pm_dis.R in Guan_and_Garrido-Oter_2025 repo
## source functions for permutation
source("functions/fun_get_spectra_and_get_dis_df.R")
source("functions/fun_get_egv.R")
source("functions/fun_get_pm_grp.R")
algae <- bs_pm(algae, group = "Cluster", 
               g_size = 30, s_size = 15, bs = 11, pm = 11)
algae <- net_dis(algae, group ="Cluster", method = "spectra")
dis_stat(algae)

out_dir <- "./Order/"
out_dir1 <- paste0(out_dir, "/single_grp_pm/")
out_dir2 <- paste0(out_dir, "/single_grp_pm_egv/")
out_dir3 <- paste0(out_dir, "/spectra_dis/")
if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(out_dir1)) dir.create(out_dir1)
if (!dir.exists(out_dir2)) dir.create(out_dir2)
if (!dir.exists(out_dir3)) dir.create(out_dir3)

grp <- tax[tax$FeatureID %in% rownames(norm(algae)), ]
grp <- grp[!(is.na(grp$Order)), ]
group <- unique(as.character(grp$Order))


len_f <- length(group)

for (i in 1:len_f) {
  grp_lst <- f_name <- group[i]
  this_grp <- grp[grp$Order %in% grp_lst, ]
  
  asv_lst <- this_grp$FeatureID
  
  this_out_dir <- paste0(out_dir1, f_name, "/")
  if (!dir.exists(this_out_dir)) dir.create(this_out_dir)
  
  print(asv_lst)
  this_pm <- get_pm_group(norm(algae), des(algae), asv_lst = asv_lst,
                          group = "Cluster", g_size = 30, s_size = 15,
                          pm = 11,
                          dir = this_out_dir)
  
  this_out_dir_egv <- paste0(out_dir2, f_name, "/")
  if (!dir.exists(this_out_dir_egv)) dir.create(this_out_dir_egv)
  
  get_egv(this_out_dir, out_dir = this_out_dir_egv)
}

################################################################################
## modified from get_dis_and_mean.R in Guan_and_Garrido-Oter_2025 repo
cmp_lst <- c("1_vs_2", "1_vs_3", "2_vs_3")
egv_folder <- "./Order/single_grp_pm_egv/"

family <- list.dirs(path = egv_folder, full.names = FALSE, recursive = FALSE)
len_f <- length(family)
dis_bs_algae <- dis_bs(algae)

for (cmp in cmp_lst){
  ### get this_bs_dis
  cmp_g1 <- unlist(strsplit(cmp, "_vs_"))[1]
  cmp_g2 <- unlist(strsplit(cmp, "_vs_"))[2]
  this_bs_dis <- dis_bs_algae[(dis_bs_algae$Group1 == cmp_g1 & dis_bs_algae$Group2 == cmp_g2)|
          (dis_bs_algae$Group1 == cmp_g2 & dis_bs_algae$Group2 == cmp_g1), ]
  
  this_cmp_test <- c()
  
  for (i in 1:len_f){
    this_dis <- paste0(egv_folder, "/", family[i], "/dis_spectra_pm_",
                       cmp, ".txt")
    #print(this_dis) 
    #next
    this_dis <- read.table(this_dis, header = T, sep = "\t")
    
    ## test if two distance are significantly different
    this_test <- t.test(this_bs_dis$Distance,
                        this_dis$Distance, alternative = "greater")
    
    #        if (this_test$p.value > 0.05) next
    m1 <- mean(this_bs_dis$Distance)
    m2 <- mean(this_dis$Distance)
    m12 <- m1 - m2
    
    this_cmp_test <- rbind(this_cmp_test,
                           c(family[i], this_test$p.value, m1, m2, m12))
  }
  
  colnames(this_cmp_test) <-c("Order", "P_value",
                              "BS_Mean", "PM_fam_Mean", "Decrease")
  write.table(this_cmp_test, paste0(out_dir3, "/", cmp, ".txt"),
              quote = F, sep = "\t", row.names = F)
}

dis_files <- list.files(out_dir3, pattern = ".txt$", full.names = TRUE)

all_dis <- c()

for (d in dis_files) {
  this_cmp <- unlist(strsplit(d, "/"))[length(unlist(strsplit(d, "/")))]
  this_cmp <- unlist(strsplit(this_cmp, ".txt"))[1]
  
  g1 <- unlist(strsplit(this_cmp, "_vs_"))[1]
  g2 <- unlist(strsplit(this_cmp, "_vs_"))[2]
  
  if (g1 == g2) next
  
  this_dis <- read.table(d, header = T, sep = "\t")
  this_dis <- this_dis[, c("Order", "Decrease", "P_value")]
  
  this_dis$Comparison <- this_cmp
  
  all_dis <- rbind(this_dis, all_dis)
}

write.table(all_dis, "Order_single_pm_dis.txt",
            quote = F, sep = "\t", row.names = F)

all_dis <- all_dis[, c("Order", "Decrease")]

all_dis_mean <- all_dis %>% group_by(Order) %>%
  summarise(Decrease_Mean = mean(Decrease))

all_dis_mean <- as.data.frame(all_dis_mean)

write.table(all_dis_mean, "Order_single_pm_dis_mean.txt",
            quote = F, sep = "\t", row.names = F)

