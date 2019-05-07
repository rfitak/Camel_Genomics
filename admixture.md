# Analysis of Genetic Structure

_LD Prune the SNP dataset in R using SNPRelate v1.10.1_
```R
# Load SNPRelate and other dependencies
source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")
library(gdsfmt)
library(SNPRelate)

# Convert VCF to GDS format
vcf.fn="All.SNPs.filtered.vcf"
snpgdsVCF2GDS(vcf.fn, "All.camels.gds", method = "biallelic.only")
genofile = snpgdsOpen("All.camels.gds")

# Get sample (camel) IDs
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# LD Prune at 0.5 level
snp.list = snpgdsLDpruning(genofile,
   sample.id = samp.id,
   autosome.only = FALSE,
   remove.monosnp = TRUE,
   maf = NaN,
   missing.rate = NaN,
   method = "corr",
   slide.max.bp = 1000000,
   slide.max.n = NA,
   ld.threshold = 0.5,
   num.thread = 16,
   verbose = TRUE)
snpset.id <- unlist(snp.list)
# Results: 90,918 SNPs selected
```

_Construct and plot Principal Components Analysis (PCA)_
```R
# Load GDS file (if necessary)
genofile = snpgdsOpen("All.camels.gds")

# Make PCA using SNPRelate v1.10.1
# Filter for the 90,918 unlinked SNPs
pca = snpgdsPCA(genofile,
   snp.id = snp.list,
   sample.id = samp.id,
   num.thread = 1,
   verbose = T,
   autosome.only = FALSE)

# Plot the results
col = c(rep("#d40000", 7),rep("#000000",9), rep("#3771c8",9))
pdf("PCA.pdf")
plot(pca$eigenvect[,2], pca$eigenvect[,1],
   xlab = "Principle Component 2",
   ylab = "Principle Component 1",
   col = col,
   pch = 19,
   cex = 1.2,
   cex.lab = 1.3,
   las = 1)
leg = c(expression(italic("C. dromedarius")),
   expression(italic("C. bactrianus")),
   expression(italic("C. ferus")))
legend("topright", legend = leg,
   col = c("#000000", "#d40000", "#3771c8"), pch = 19, bty = "y")
# text(pca$eigenvect[,2], pca$eigenvect[,1], labels = samp.id)        # check sample IDs
dev.off()

pdf("PCA-bars.pdf")
pc.percent <- pca$varprop*100
barplot(pc.percent[1:10], axes = F, ylim = c(0, 25), names.arg = c(1:10), col = "grey35")
axis(2, las = 1, at = c(0, 5, 10, 15, 20, 25), labels = c(0, 5, 10, 15, 20, 25))
title(xlab = "Principle Component", ylab = "Percent Variation", cex.lab = 1.3)
dev.off()

pdf("PCA-pairs.pdf")
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits = 2), "%", sep = "")
pairs(pca$eigenvect[,1:4], col = col, labels = lbls, pch = 19, las = 1, cex.axis = 1.3)
dev.off()
```
