# Window-based analysis of genetic variation.
The following code is for calculating various measures of genetic variation across sliding windows.  The code is based on python scripts to manipulate window-based approaches on VCF files of genomic variant data.  The python scripts are hosted on [Simon Martin's GitHub Page](https://github.com/simonhmartin/genomics_general).

First, install Simon Martin's 'genomics_general' toolset:
```bash
git clone https://github.com/simonhmartin/genomics_general.git

# Note: I am using OSX v10.13.6; python v2.7.15; numpy v1.15.4
```

Next, we need to prepare our VCF file of SNPs by compression and index it.  This should be done using the [Samtools](http://www.htslib.org) suite of tools.
```bash
bgzip -c snps.vcf > snps.vcf.gz
tabix -p vcf snps.vcf.gz
```

Next, we convert the VCF formatted genotypes to the .geno format required by the python scripts.  The format is quite simple, and can be modified as needed.  See the manual for the python scripts by running them with the ```-h``` option.  For example:
```bash
python genomics_general/VCF_processing/parseVCF.py -h
```
Convert the VCF file.
```bash
python genomics_general/VCF_processing/parseVCF.py \
   --infile snps.vcf.gz \
   --outfile input.geno.gz \
   --excludeFile XY.exclude
```
Notes:
- --exclude XY.exclude :: a file listing the scaffold IDs to exclude.  These are putative X and Y scaffolds identified previously - see file here [XY.exclude](./Data/XY.exclude)
-  If you add ".gz" to the output file, it is automatically compressed using gzip.
- We began with 10819573 SNPs, but reduced to 10581849 SNPs after removing the X and Y scaffold loci

Now we are ready to calculate the various metrics!!!
The following script will be repeated multiple times using the ```--analysis``` flag to calculate various measures of genetic variation and differentiation.  We will combine all the data at the end into a single file of results.
- *__popFreq__* will calculate nucleotide diversity (pi), Watterson's theta, and Tajima's D for each population
- *__indHet__* will calculate heterozygosity for each individual
- *__popDist*__* will calculate xxxx for each individual
- *__popPairDist__* will calculate divergence (Dxy) and (Fst) for each population pair

```bash
# Must be run from the same location as the genomics.py module

# popFreq
python2.7 popgenWindows.py \
   --windType coordinate \
   --windSize 100000 \
   --stepSize 50000 \
   --minSites 10 \
   --writeFailedWindows \
   --genoFile ../input.geno.gz \
   --analysis popFreq \
   --roundTo 4 \
   --outFile popFreq.100kb.csv.gz \
   --genoFormat phased \
   --exclude ../XY.exclude \
   --Threads 2 \
   --verbose \
   --population Drom DROM802,Drom439,Drom795,Drom796,Drom797,Drom800_55,Drom806,Drom816,Drom820 \
   --population DC DC158_53,DC269,DC399,DC400,DC402,DC408,DC423 \
   --population WC WC214,WC216,WC218,WC219,WC220,WC247,WC303_108,WC304,WC305

# indHet
python2.7 popgenWindows.py \
   --windType coordinate \
   --windSize 100000 \
   --stepSize 50000 \
   --minSites 10 \
   --writeFailedWindows \
   --genoFile ../input.geno.gz \
   --analysis indHet \
   --roundTo 4 \
   --outFile indHet.100kb.csv.gz \
   --genoFormat phased \
   --exclude ../XY.exclude \
   --Threads 2 \
   --verbose \
   --population Drom DROM802,Drom439,Drom795,Drom796,Drom797,Drom800_55,Drom806,Drom816,Drom820 \
   --population DC DC158_53,DC269,DC399,DC400,DC402,DC408,DC423 \
   --population WC WC214,WC216,WC218,WC219,WC220,WC247,WC303_108,WC304,WC305

# popDist
python2.7 popgenWindows.py \
   --windType coordinate \
   --windSize 100000 \
   --stepSize 50000 \
   --minSites 10 \
   --writeFailedWindows \
   --genoFile ../input.geno.gz \
   --analysis popDist \
   --roundTo 4 \
   --outFile popDist.100kb.csv.gz \
   --genoFormat phased \
   --exclude ../XY.exclude \
   --Threads 2 \
   --verbose \
   --population Drom DROM802,Drom439,Drom795,Drom796,Drom797,Drom800_55,Drom806,Drom816,Drom820 \
   --population DC DC158_53,DC269,DC399,DC400,DC402,DC408,DC423 \
   --population WC WC214,WC216,WC218,WC219,WC220,WC247,WC303_108,WC304,WC305

# popPairDist
python2.7 popgenWindows.py \
   --windType coordinate \
   --windSize 100000 \
   --stepSize 50000 \
   --minSites 10 \
   --writeFailedWindows \
   --genoFile ../input.geno.gz \
   --analysis popPairDist \
   --roundTo 4 \
   --outFile popPairDist.100kb.csv.gz \
   --genoFormat phased \
   --exclude ../XY.exclude \
   --Threads 2 \
   --verbose \
   --population Drom DROM802,Drom439,Drom795,Drom796,Drom797,Drom800_55,Drom806,Drom816,Drom820 \
   --population DC DC158_53,DC269,DC399,DC400,DC402,DC408,DC423 \
   --population WC WC214,WC216,WC218,WC219,WC220,WC247,WC303_108,WC304,WC305
   
# Old version...
python genomics_general/popgenWindows.py \
   --windType coordinate \
   --windSize 100000 \
   --stepSize 50000 \
   --minSites 10 \
   --genoFile input.geno.gz \
   --outFile variation.100kb.csv.gz \
   --genoFormat phased \
   --indHet \
   --exclude xy.exclude \
   --Threads 8 \
   --verbose \
   --population Drom DROM802,Drom439,Drom795,Drom796,Drom797,Drom800_55,Drom806,Drom816,Drom820 \
   --population DC DC158_53,DC269,DC399,DC400,DC402,DC408,DC423 \
   --population WC WC214,WC216,WC218,WC219,WC220,WC247,WC303_108,WC304,WC305
```
Parameters explained:
- --windType coordinate :: a coordinate-based window length rather than number of SNPs
- --windSize 100000 :: 100kb window length
- --stepSize 50000 :: window moves in 50kb steps
- --minSites 10 :: at least 10 SNPs must be in a window to calcualte statistics
- --genoFormat phased :: input format, only change if you set this in the previous step to be different
- --writeFailedWindows :: record windows that fail for having less than 10 sites
- --analysis popFreq :: choise of which measurements to record
- --roundTo 4 :: round to four digits
- --exclude XY.exclude :: remove scaffolds in this file (although it was also performed in the previous step)
- --Threads 2 :: use 2 threads, this took only ~15 minutes using 2 threads
- --verbose :: provide detailed run time log
- --population XX 1,2,3,... :: the population ID and comma-separated list of samples.


Last, we will combine these results together in a final data table:
```bash
paste \
   -d"," \
   <(gunzip -c popFreq.100kb.csv.gz)
   <gunzip -c indHet.100kb.csv.gz | cut -d"," -f6-) \
   <gunzip -c popDist.100kb.csv.gz | cut -d"," -f6-) \
   <gunzip -c popPairDist.100kb.csv.gz | cut -d"," -f6-) | \
   gzip > Final.100kb.csv.gz
```


### Find windows in Drom and WC under positive selection
Here we will attempt to find windows with a dearth or polymorphism and excess of divergence with wild camels.  We will use a cutoff of 99.5% and 0.5% for the extreme 'outliers'.  Later, we will add the results of the population branch statistic to make one big file!
```R
# First load the window results file
input = read.csv(gzfile("popFReq.100kb.csv.gz"), header = T)
   # 39660 windows in total
   
# Remove windows with less than 10 sites
data = subset(input, sites >= 10)
   # 37030 windows passed (removed 2630 windows)

# Find Drom outlier windows
Drom.selected = data[which(data$thetaPi_Drom < quantile(data$thetaPi_Drom, 0.005) & data$dxy_Drom_WC > quantile(data$dxy_Drom_WC, 0.995)), ]
write.table(Drom.selected, file = "Drom.selected.100kb.tsv", quote = F, sep = "\t", row.names = F)

# Find DC outlier windows
DC.selected = data[which(data$pi_DC < quantile(data$pi_DC, 0.005) & data$dxy_DC_WC > quantile(data$dxy_DC_WC, 0.995)), ]
write.table(DC.selected, file = "DC.selected.100kb.tsv", quote = F, sep = "\t", row.names = F)


```





### Find windows under positive selection
In this section we will find the extreme outlier windows to identify regions under positive selection.  We will also add the calculation of the population branch statistic (PBS) to further identify positive selection in the DC or WC lineages.  See [Yi et al. 2010 Science](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3711608/) for the original derivation.  It is actually quite simple.
Most of this analysis will be done in R.

The function below calculates the PBS statistic.  The input is a vector of three Fst measures. The first Fst is the population of interest compared with the sister population, the second is population of interest compared with the outgroup, and the last Fst is the sister population compared with the outgroup.  The function returns the PBS for the population of interest.

```R
PBS = function(fst){
   if (any(is.na(fst) == T)) return(NA)
   else {
      if (fst[1] == 1) fst[1] = 0.9999
      if (fst[2] == 1) fst[2] = 0.9999
      if (fst[3] == 1) fst[3] = 0.9999
      T1 = -log(1 - fst[1]) #DC vs WC
      T2 = -log(1 - fst[2]) #DC vs DROM
      T3 = -log(1 - fst[3]) #WC cs DROM
      pbs = (T1 + T2 - T3) / 2
      return(pbs)
   }
}
```

### Calculate Reynold's Fst
Unfortunately, the PBS is defined in terms of an Fst value calculated according to the [Reynolds et al (1983)](http://image.sciencenet.cn/olddata/kexue.com.cn/blog/admin/images/upfiles/2008423153436872569.pdf) method. The Fst calculated by the python tool above is slightly different.  Here is the R code to do so:

```R
# First load the genotype file and window results file
genos = read.table(gzfile("input.geno.gz"), sep = "\t", header = T, comment.char = "", stringsAsFactors = F)
data = read.csv(gzfile("variation.100kb.csv.gz"), header = T)

# Make list of individuals in each population
Drom = c("DROM802", "Drom439", "Drom795", "Drom796", "Drom797", "Drom800_55", "Drom806", "Drom816", "Drom820")
DC = c("DC158_53", "DC269", "DC399", "DC400", "DC402", "DC408", "DC423")
WC = c("WC214", "WC216", "WC218", "WC219", "WC220", "WC247", "WC303_108", "WC304", "WC305")

# Setup empty vector of new Reynolds' Fst
Reynolds.fst = vector()

# Now we will cycle through each window in 'data'
for (w in 1:nrow(data)){
   # Get snps in the first window
   chrom = as.character(data$scaffold[w])
   start = data$start[w]
   end = data$end[w]
   window = genos[which(genos$X.CHROM == chrom),]
   window = window[window$POS >= start & window$POS <= end,]

   # Get major and minor allele overall for each locus
   window.tmp = window[,3:ncol(window)]
   alleles = apply(window.tmp, 1, function(x) gsub("/", "", paste(x[which(x != "N/N")], collapse="")))
   counts = lapply(strsplit(alleles, ""), function(x) summary(as.factor(x)))

   # Subsets the window data by population
   Drom.win = window[,which(colnames(window) %in% Drom)]
   DC.win = window[,which(colnames(window) %in% DC)]
   WC.win = window[,which(colnames(window) %in% WC)]
   Drom.alleles = strsplit(apply(Drom.win, 1, function(x) gsub("/", "", paste(x[which(x != "N/N")], collapse=""))), "")
   DC.alleles = strsplit(apply(DC.win, 1, function(x) gsub("/", "", paste(x[which(x != "N/N")], collapse=""))), "")
   WC.alleles = strsplit(apply(WC.win, 1, function(x) gsub("/", "", paste(x[which(x != "N/N")], collapse=""))), "")

   # Now process each SNP in the window to get Reynolds' a and a+b
   out = vector()
   for (i in 1:nrow(window)){
      Drom.counts = sapply(names(counts[[i]]), function(x) length(which(Drom.alleles[[i]] == x)))
      Drom.freqs = Drom.counts / sum(Drom.counts)
      Drom.N = sum(Drom.counts)
      Drom.alpha1 = 1 - sum((Drom.freqs)^2)

      DC.counts = sapply(names(counts[[i]]), function(x) length(which(DC.alleles[[i]] == x)))
      DC.freqs = DC.counts / sum(DC.counts)
      DC.alpha1 = 1 - sum((DC.freqs)^2)
      DC.N = sum(DC.counts)

      WC.counts = sapply(names(counts[[i]]), function(x) length(which(WC.alleles[[i]] == x)))
      WC.freqs = WC.counts / sum(WC.counts)
      WC.alpha1 = 1 - sum((WC.freqs)^2)
      WC.N = sum(WC.counts)

      # For each pair of populations
      # Drom v DC
      Drom.vs.DC.al = 1/2 * sum((Drom.freqs - DC.freqs)^2) - (Drom.N + DC.N) * (Drom.N * Drom.alpha1 + DC.N * DC.alpha1) / (4 * Drom.N * DC.N * (Drom.N + DC.N - 1))
      Drom.vs.DC.abl = 1/2 * sum((Drom.freqs - DC.freqs)^2) + (4 * Drom.N * DC.N - Drom.N - DC.N) * (Drom.N * Drom.alpha1 + DC.N * DC.alpha1) / (4 * Drom.N * DC.N * (Drom.N + DC.N - 1))
   
      # Drom v WC
      Drom.vs.WC.al = 1/2 * sum((Drom.freqs - WC.freqs)^2) - (Drom.N + WC.N) * (Drom.N * Drom.alpha1 + WC.N * WC.alpha1) / (4 * Drom.N * WC.N * (Drom.N + WC.N - 1))
      Drom.vs.WC.abl = 1/2 * sum((Drom.freqs - WC.freqs)^2) + (4 * Drom.N * WC.N - Drom.N - WC.N) * (Drom.N * Drom.alpha1 + WC.N * WC.alpha1) / (4 * Drom.N * WC.N * (Drom.N + WC.N - 1))
   
      # DC v WC
      DC.vs.WC.al = 1/2 * sum((DC.freqs - WC.freqs)^2) - (DC.N + WC.N) * (DC.N * DC.alpha1 + WC.N * WC.alpha1) / (4 * DC.N * WC.N * (DC.N + WC.N - 1))
      DC.vs.WC.abl = 1/2 * sum((DC.freqs - WC.freqs)^2) + (4 * DC.N * WC.N - DC.N - WC.N) * (DC.N * DC.alpha1 + WC.N * WC.alpha1) / (4 * DC.N * WC.N * (DC.N + WC.N - 1))

      # Add to output vector
      out = rbind(out, c(Drom.vs.DC.al, Drom.vs.DC.abl, Drom.vs.WC.al, Drom.vs.WC.abl, DC.vs.WC.al, DC.vs.WC.abl))
   }
   
   # Now get weighted Fst per locus, or ratio of the averages of al and abl
   Drom.vs.DC.fst = mean(out[,1], na.rm = T) / mean(out[,2], na.rm = T)
   Drom.vs.WC.fst = mean(out[,3], na.rm = T) / mean(out[,4], na.rm = T)
   DC.vs.WC.fst = mean(out[,5], na.rm = T) / mean(out[,6], na.rm = T)
   Reynolds.fst = rbind(Reynolds.fst, c(Drom.vs.DC.fst, Drom.vs.WC.fst, DC.vs.WC.fst))
   
   # Calculate progress of iterations
   p = w / nrow(data)
   if (p %in% seq(0,1, by = 0.1)) message(paste0(p*100, "% complete..."))
}
```

So now we have a data frame ```Reynolds.fst``` of Reynolds' Fst values for each population pair in each window.  Finally, we will calculate the PBS for both DC and WC and append our results to the end of our data file.
```R
# Add column names to Reynolds.fst and merge with data
colnames(Reynolds.fst) = c("Rey_Fst_Drom_DC", "Rey_Fst_Drom_WC", "Rey_Fst_DC_WC")
data = cbind(data, Reynolds.fst)

# Get the PBS for DC
fst = data[,c(42, 40, 41)]
PBS_DC = apply(fst, 1, PBS)

# Get the PBS for WC
fst = data[,c(42, 41, 40)]
PBS_WC = apply(fst, 1, PBS)

# Merge with data
data = cbind(data, PBS_DC, PBS_WC)

# Save the new data table
write.table(data, file = "Final-table.100kb.tsv", quote = F, sep = "\t", row.names = F)

# Find Drom outlier windows
Drom.selected = data[which(data$pi_Drom < quantile(data$pi_Drom, 0.005) & data$dxy_Drom_WC > quantile(data$dxy_Drom_WC, 0.995)), ]
write.table(Drom.selected, file = "Drom.selected.100kb.tsv", quote = F, sep = "\t", row.names = F)

# Find DC outlier windows
DC.selected = data[which(data$pi_DC < quantile(data$pi_DC, 0.005) & data$dxy_DC_WC > quantile(data$dxy_DC_WC, 0.995)), ]
write.table(DC.selected, file = "DC.selected.100kb.tsv", quote = F, sep = "\t", row.names = F)

# DC PBS outlier windows
DC.PBS = data[which(data$PBS_DC > quantile(data$PBS_DC, 0.999, na.rm = T)), ]
write.table(DC.PBS, file = "DC.PBS.100kb.tsv", quote = F, sep = "\t", row.names = F)

# WC PBS outlier windows
WC.PBS = data[which(data$PBS_WC > quantile(data$PBS_WC, 0.999, na.rm = T)), ]
write.table(WC.PBS, file = "WC.PBS.100kb.tsv", quote = F, sep = "\t", row.names = F)
```

YES!!!!!!  We now have files containing the regions under selection!!!!
The following commands can find the genes in these regions.
```bash
# Drom
bedtools \
   intersect \
   -wa \
   -a /wrk/rfitak/DONOTREMOVE/SNP-ANALYSIS/ref_CB1_scaffolds.gff3 \
   -b <(cut -f1-3 Drom.selected.100kb.tsv | sed '1d') | \
   grep -P "\tgene\t" | \
   sort | \
   uniq > Drom.selected.genes.gff3
   # Results: 5 tRNAs, 3 pseudogenes, 23 genes!

# DC
bedtools \
   intersect \
   -wa \
   -a /wrk/rfitak/DONOTREMOVE/SNP-ANALYSIS/ref_CB1_scaffolds.gff3 \
   -b <(cut -f1-3 DC.selected.100kb.tsv | sed '1d') | \
   grep -P "\tgene\t" | \
   sort | \
   uniq > DC.selected.genes.gff3
   # Results: no genes overlapping

# DC PBS
bedtools \
   intersect \
   -wa \
   -a /wrk/rfitak/DONOTREMOVE/SNP-ANALYSIS/ref_CB1_scaffolds.gff3 \
   -b <(cut -f1-3 DC.PBS.100kb.tsv | sed '1d') | \
   grep -P "\tgene\t" | \
   sort | uniq > DC.PBS.genes.gff3
   # Results: 2 tRNAs, 10 genes!

# WC PBS
bedtools \
   intersect \
   -wa \
   -a /wrk/rfitak/DONOTREMOVE/SNP-ANALYSIS/ref_CB1_scaffolds.gff3 \
   -b <(cut -f1-3 WC.PBS.100kb.tsv | sed '1d') | \
   grep -P "\tgene\t" | \
   sort | uniq > WC.PBS.genes.gff3
   # Results: 19 tRNAs, 1 pseudogenes, 64 genes!
```



# Old Code
Sliding window-based anlaysis using the [PopGenome](https://cran.r-project.org/web/packages/PopGenome/index.html) package in R
```R
library(PopGenome)
filename = "snps.vcf.gz"
vcf_handle <- .Call("VCF_open",filename)
samplenames <- .Call("VCF_getSampleNames",vcf_handle)

GENOME.class <- readVCF("snps.vcf.gz",numcols = 100000, tid = "NW_006223456.1", 
   from = 1, to = 1000000000, include.unknown = F, approx = FALSE, out = "", 
   parallel = FALSE, gffpath = FALSE)
DC = samplenames[1:7]
DROM = samplenames[8:16]
WC = samplenames[17:25]
GENOME.class <- set.populations(GENOME.class,list(DC, DROM, WC), diploid = TRUE)
#GENOME.class <- set.outgroup(GENOME.class,c("z"), diploid = TRUE)

slide <- sliding.window.transform(GENOME.class,10000,10000, type=2)
length(slide@region.names)
slide <- diversity.stats(slide)
```
Sliding window-based anlaysis using [genomics.py](https://github.com/simonhmartin/genomics_general/blob/master/genomics.py) from Simon Martin
```python
import gzip
import genomics
genoFileName = "input2.geno.gz"
genoFile = gzip.open(genoFileName, "r")
popNames = ["Cdrom","Cbact","Cferu"]
samples = [["Cdrom.DROM802","Cdrom.Drom439","Cdrom.Drom795","Cdrom.Drom796",
           "Cdrom.Drom797","Cdrom.Drom800_55","Cdrom.Drom806","Cdrom.Drom816","Cdrom.Drom820"],
           ["Cbact.DC158_53","Cbact.DC269","Cbact.DC399","Cbact.DC400","Cbact.DC402","Cbact.DC408","Cbact.DC423"],
           ["Cferu.WC214","Cferu.WC216","Cferu.WC218","Cferu.WC219","Cferu.WC220",
           "Cferu.WC247","Cferu.WC303_108","Cferu.WC304","Cferu.WC305"]]
samples = [["DROM802","Drom439","Drom795","Drom796",
           "Drom797","Drom800_55","Drom806","Drom816","Drom820"],
           ["DC158_53","DC269","DC399","DC400","DC402","DC408","DC423"],
           ["WC214","WC216","WC218","WC219","WC220",
           "WC247","WC303_108","WC304","WC305"]]

sampleData = genomics.SampleData(popInds = samples, popNames = popNames)

windSize = 10000
stepSize = 5000
minSites = 2

windowGenerator = genomics.slidingWindows(genoFile, windSize, stepSize, skipDeepcopy = True)
outFileName = "pop_div.w50m25s25.csv"
outFile = open(outFileName, "w")
outFile.write("scaffold,start,end,sites")

for x in range(len(popNames)-1):
    outFile.write(",pi_" + popNames[x])
    for y in range(x+1,len(popNames)):
        outFile.write(",Fst_" + popNames[x] + "_" + popNames[y])
        outFile.write(",dxy_" + popNames[x] + "_" + popNames[y])

outFile.write("\n")

n=0
for window in windowGenerator:
    if window.seqLen() >= minSites:
        #if there are enough sites, make alignment object
        Aln = genomics.callsToAlignment(window.seqDict(), sampleData, seqType = "pairs")
        #get divergence stats
        statsDict = genomics.popDiv(Aln)
        stats = []
        for x in range(len(popNames)-1):
            #retrieve pi for each population
            stats.append(statsDict["pi_" + popNames[x]])
            for y in range(x+1,len(popNames)):
                #retrieve dxy and Fst for each pair
                stats.append(statsDict["fst_" + popNames[x] + "_" + popNames[y]])
                stats.append(statsDict["dxy_" + popNames[x] + "_" + popNames[y]])
        #add window stats and write to file
        out_data = [window.scaffold, window.start, window.end, window.seqLen()] + [round(s,3) for s in stats]
        outFile.write(",".join([str(x) for x in out_data]) + "\n")
    print n
    n+=1
outFile.close()
genoFile.close()
exit()
```
