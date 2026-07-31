##rarefaction curves datasets
asv_algae<- read.table("otu_table.txt", header=TRUE, sep="\t")

asv_bac<- read.table("asv_table_bac.txt", header=TRUE, sep="\t")

asv_fun<- read.table("asv_table_fun.txt", header=TRUE, sep="\t")

library(vegan)

# rows = samples, cols = taxa
rarecurve(asv_algae, step = 500, col = "black", label = FALSE, xlim = c(0, 10000))
abline(v = 2000, col = "red", lty = 2, lwd = 2)

rarecurve(asv_bac, step = 500, col = "black", label = FALSE, xlim = c(0, 10000))
abline(v = 2000, col = "red", lty = 2, lwd = 2)

rarecurve(asv_fun, step = 500, col = "black", label = FALSE, xlim = c(0, 10000))
abline(v = 2000, col = "red", lty = 2, lwd = 2)
