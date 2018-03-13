# GO Term Enrichment using TopGO
This section describes the basic code for producing GO term enrichment analyses in R using the Bioconductor package [topGO](http://bioconductor.org/packages/release/bioc/html/topGO.html). To run TopGO, you will need:
  1.  A GO term mapping database the lists each GO term assigned to each gene
  2.  A list of reference genes
  3.  A list of test genes

First, download and install topGO in R
```R
# Download and Install topGO
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("topGO")
library(topGO)
```

Next, we made a GO mapping file from the Blast2GO annotations in the format:  
GeneID   GO, GO, GO  
This mapping file can be found here: [gene2GO](./Data/gene2GO).  
There are 15,603 genes with GO annotations, and 17,175 genes in total in the reference set.  The reference list of proteins can be found here: [reference-set.protein.IDs](./Data/reference-set.protein.IDs).

In R, we now load our reference gene list, Go information, and build a function to pick selected genes (simply 1 or 0 for selected or not, respectively).
```R
# Load gene universe (only the 4th column to jsut grab the refseq ID)
genes.universe = as.vector(read.table("reference-set.protein.IDs", header = F, sep = "|")[,4])

# Build selection function
picker=function (allScore) {
    return(allScore == 1)
}

# Load GO database mappings
geneID2GO <- readMappings(file = "gene2GO", sep="\t", IDsep=",")

# Note: to change GO mappings to a list of GO terms with genes mapping to each:
GO2geneID <- inverseList(geneID2GO)
```

Great!  We have the reference data loaded and ready for use in topGO.  Next, we setup the information for the type of test we want to perform (here we just use the classic Fisher exact test).  The 'weight01' test is also recommended but beware issues with p-values since all tests are not independent.  We also load our set of test protein IDs (this file varies between enrichment analyses and only an example is shown).
```R
# Setup test statistic to output
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

# Load list of selected genes
genes.selected = scan("selected.genes", what = "character")

# Label all genes as 0/1 for nonselected/selected
geneList = as.numeric(as.character(factor(as.integer (genes.universe %in% genes.selected))))

# Assign names to above
names(geneList) <- genes.universe
```
Finally, we build and run our topGO dataset and generate a final results table with FDR corrected p-values.  
```R
# Make topGO variable (repeat for "BP", "MF", or "CC" GO types)
GOdata <- new("topGOdata",
description = "Selection Run 1",
ontology = "BP",
allGenes = geneList,
geneSelectionFun = picker,
nodeSize = 5,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

# Load all GO terms used to adjust later
allGO = usedGO(object = GOdata)

# Run enrichment test
resultFisher <- getSigGroups(GOdata, test.stat)
   
# Create table
sig.tab <- GenTable(GOdata, Fis = resultFisher, topNodes = length(allGO))

# Adjust p-values
AdjP=p.adjust(sig.tab$Fis, method = "fdr")
   sig.tab=cbind(sig.tab, AdjP)
   
# Write output
write.table(sig.tab, file = "Selection_Run1.tsv", sep = "\t", quote = F, row.names = F)
```
