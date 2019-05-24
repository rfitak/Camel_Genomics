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

_Convert GDS file to LD-pruned, PLINK-formatted bed|bim|fam files for use in Admixture_
```R
# Load libraries (if not already loaded
library(gdsfmt)
library(SNPRelate)

# First make a GDS file as a subset of the first (LD-pruned GDS file)
snpgdsCreateGenoSet(src.fn = "All.camels.gds",
   dest.fn = "test.gds",
   sample.id = samp.id,
   snp.id = snp.list[,1],
   snpfirstdim = NULL,
   compress.annotation = "ZIP.max",
   compress.geno = "",
   verbose = TRUE)

# Then open the new GDS file
genofile = snpgdsOpen("test.gds")

# Convert chromosome (scaffold) names to numbers
x = as.numeric(as.factor(read.gdsn(index.gdsn(genofile, "snp.chromosome"))))

# Create new GDS file
geno = read.gdsn(index.gdsn(genofile, "genotype"))
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
snp.id = read.gdsn(index.gdsn(genofile, "snp.id"))
snp.position = read.gdsn(index.gdsn(genofile, "snp.position"))
snp.allele = read.gdsn(index.gdsn(genofile, "snp.allele"))
snpgdsCreateGeno("test2.gds",
   genmat = t(geno),
   sample.id = sample.id,
   snp.id = snp.id,
   snp.chromosome = x,
	snp.position = snp.position,
   snp.allele = snp.allele)

# Open new GDS file
genofile = snpgdsOpen("test2.gds")

# Convert to PLINK format
snpgdsGDS2BED(genofile,
   bed.fn = "admixture.input",
   sample.id = samp.id,
   snpfirstdim = T,
   verbose = TRUE)
```

_Run Admixture_
```bash
# Run ADMIXTURE for K=1-5
for i in {1..5}
do
admixture \
   -j8 \
   --cv=10 \
   -B100 \
   -C 0.0001 \
   -c 0.0001 \
   admixture.input.bed $i | \
   tee log.$i.out
echo "Finished $i"
done
```

_Plot Admixture results (Cross-validation plot and ancestry plots)_
```R
# Ancestry barplots for K=2 and 3
ID = c("DC158", "DC269", "DC399", "DC400", "DC402", "DC408", "DC423",
   "WC214", "WC216", "WC218", "WC219", "WC220", "WC247", "WC303", "WC304", "WC305",
   "Drom802", "Drom439", "Drom795", "Drom796", "Drom797", "Drom800", "Drom806", "Drom816", "Drom820")
tbl2 = as.matrix(read.table("admixture.input.2.Q"))
tbl2 = rbind(tbl2[1:7,], tbl2[17:25,], tbl2[8:16,])
tbl3 = as.matrix(read.table("admixture.input.3.Q"))
tbl3 = rbind(tbl3[1:7,], tbl3[17:25,], tbl3[8:16,])

# Plot bars
pdf("Admixture.pdf")
par(mfrow = c(2, 1))
barplot(t(tbl2),
   col = c("#000000", "#d40000"),
   ylab = "Ancestry",
   border = NA,
   las = 2,
   cex.lab = 1.3)
barplot(t(tbl3),
   col = c("#d40000", "#3771c8", "#000000"),
   ylab = "Ancestry",
   border = NA,
   names.arg = ID,
   las = 2,
   cex.lab = 1.3)
dev.off()

# Plot CV error
pdf("CV_error.pdf")
x = c(1:5)
y = c(0.16216, 0.12124, 0.12273, 0.14973, 0.18435) # CV error values
plot(x, y,
   type = "l",
   lwd = 2,
   las = 1,
   xlab = expression(italic("k")),
   ylab = "CV Error",
   cex.lab = 1.3)
dev.off()
```
