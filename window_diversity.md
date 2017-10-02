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
In this section we will find the extreme outlier windows to identify regions under positive selection.  Most of this analysis will be done in R.
