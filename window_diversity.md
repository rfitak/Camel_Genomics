# Window-based analysis of genetic variation.
The following code is for calculating various measures of genetic variation across sliding windows.  The code is based on python scripts to manipulate window-based approaches on VCF files of genomic variant data.  The python scripts are hosted on [Simon Martin's GitHub Page](https://github.com/simonhmartin/genomics_general).

First, install Simon Martin's 'genomics_general' toolset:
```bash
git clone https://github.com/simonhmartin/genomics_general.git
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
   --excludeFile xy.exclude
```
Notes:
- --exclude xy.exclude :: a file listing the scaffold IDs to exclude.  These are putative X and Y scaffolds identified previously
-  If you add ".gz" to the output file, it is automatically compressed using gzip.
- We began with 10819573 SNPs, but reduced to 10581849 SNPs after removing the X and Y scaffold loci

The following script will calculate individual heterozygosity, pi (within populations) and divergence (Dxy, between populations).
```bash
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
- --indHet :: calculate individual level heterozygosity
- --exclude xy.exclude :: remove scaffolds in this file (although it was also performed in the previous step)
- --Threads 8 :: use 8 threads, this took only 9 minutes using 8 threads
- --verbose :: provide detailed run time log
- --population XX 1,2,3,... :: the population ID and comma-separated list of samples.


### Find windows under positive selection
In this section we will find the extreme outlier windows to identify regions under positive selection.  We will also add the calculation of the population branch statistic (PBS) to further identify positive selection in the DC or WC lineages.  See XXXXX for the original derivation.  It is actually quite simple.
Most of this analysis will be done in R.

### Calculate Reynold's Fst
In order to calculate the PBS, we first need to calculate Fst using the equation by Reynolds et al. 1983.  The Fst calculated by the python tool above is slightly different.  Here is the R code to do so:
```R
# Reynolds Fst



```

### PBS function in R
This function calculates the PBS statistic.  The input is a vector of three Fst measures (Reynold's method). The first Fst is the population of interest compared with the sister population, the second is population of interest compared with the outgroup, and the last Fst is the sister population compared with the outgroup.  The function returns the PBS for the population of interest.

```R
PBS = function(fst){
   if (any(is.na(fst) == T)) return(NA)
   else {
      T1 = -log(1 - fst[1]) #DC vs WC
      T2 = -log(1 - fst[2]) #DC vs DROM
      T3 = -log(1 - fst[3]) #WC cs DROM
      pbs = (T1 + T2 - T3) / 2
      return(pbs)
   }
}
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
