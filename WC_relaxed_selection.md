# Code for detecting relaxed selection in the domesticated camels (dromedary and domestic Bactrian) compared with the wild camel.

### Step 1:  Calculate Tajima's D in 100 kb windows across the genome in _C. ferus_

```bash
# Download genomics windows tools
   # OSX v10.13.6; python v2.7.15; numpy v1.15.4
git clone https://github.com/simonhmartin/genomics_general.git

# Compress and index
~/bin/bin/bgzip -c All.SNPs.filtered.vcf > snps.vcf.gz
~/bin/bin/tabix -p vcf snps.vcf.gz

# Convert VCF file
python genomics_general/VCF_processing/parseVCF.py \
   --infile snps.vcf.gz \
   --outfile input.geno.gz \
   --excludeFile XY.exclude

python genomics_general/popgenWindows.py \
   --windType coordinate \
   --windSize 100000 \
   --stepSize 50000 \
   --minSites 10 \
   --writeFailedWindows \
   --genoFile input.geno.gz \
   --analysis popFreq \
   --roundTo 4 \
   --outFile popFreq.100kb.csv.gz \
   --genoFormat phased \
   --exclude XY.exclude \
   --Threads 8 \
   --verbose \
   --population Drom DROM802,Drom439,Drom795,Drom796,Drom797,Drom800_55,Drom806,Drom816,Drom820 \
   --population DC DC158_53,DC269,DC399,DC400,DC402,DC408,DC423 \
   --population WC WC214,WC216,WC218,WC219,WC220,WC247,WC303_108,WC304,WC305

# Tajima's D
~/bin/bin/vcftools \
   --gzvcf snps.vcf.gz \
   --remove-indv DC158_53 \
   --remove-indv DC269 \
   --remove-indv DC399 \
   --remove-indv DC400 \
   --remove-indv DC402 \
   --remove-indv DC408  \
   --remove-indv DC423 \
   --remove-indv DROM802 \
   --remove-indv Drom439 \
   --remove-indv Drom795 \
   --remove-indv Drom796 \
   --remove-indv Drom797 \
   --remove-indv Drom800_55 \
   --remove-indv Drom806 \
   --remove-indv Drom816 \
   --remove-indv Drom820 \
   --exclude-bed XY.exclude \
   --mac 1 \
   --TajimaD 100000 \
   --out WC.tajimaD
```

```R
# R script for calculating Tajima's D from a table of nucleotide diversity (pi) and # segragating sites



```
