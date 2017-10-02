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
- We began with 10819573 SNPs, but reduced to XXXX SNPs after removing the X and Y scaffold loci
