# SNP Calling Pipeline

_Call SNPs using the GATK v3.1-1 HaplotypeCaller_

```bash
# SNPs and variants are called across all samples simultaneously
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T HaplotypeCaller \
   -R CB1.fasta \
   -I Drom439.sorted.rmdup.mq20.real.bam
   -I Drom795.sorted.rmdup.mq20.real.bam
   -I Drom796.sorted.rmdup.mq20.real.bam
   -I Drom797.sorted.rmdup.mq20.real.bam
   -I Drom800.sorted.rmdup.mq20.real.bam
   -I Drom802.sorted.rmdup.mq20.real.bam
   -I Drom806.sorted.rmdup.mq20.real.bam
   -I Drom816.sorted.rmdup.mq20.real.bam
   -I Drom820.sorted.rmdup.mq20.real.bam
   -I DC158.sorted.rmdup.mq20.real.bam
   -I DC269.sorted.rmdup.mq20.real.bam
   -I DC399.sorted.rmdup.mq20.real.bam
   -I DC400.sorted.rmdup.mq20.real.bam
   -I DC402.sorted.rmdup.mq20.real.bam
   -I DC408.sorted.rmdup.mq20.real.bam
   -I DC423.sorted.rmdup.mq20.real.bam
   -I WC214.sorted.rmdup.mq20.real.bam
   -I WC216.sorted.rmdup.mq20.real.bam
   -I WC218.sorted.rmdup.mq20.real.bam
   -I WC219.sorted.rmdup.mq20.real.bam
   -I WC220.sorted.rmdup.mq20.real.bam
   -I WC247.sorted.rmdup.mq20.real.bam
   -I WC303.sorted.rmdup.mq20.real.bam
   -I WC304.sorted.rmdup.mq20.real.bam
   -I WC305.sorted.rmdup.mq20.real.bam
   -nct 32 \
   --genotyping_mode DISCOVERY \
   -stand_emit_conf 20 \
   -stand_call_conf 20 \
   -dcov 200 \
   -o GATKHC1.vcf
   
   # Result: 16,682,481 SNPs, 2,823,343 INDELS, and 19,505,824 total
```
<div style="page-break-after: always;"></div>

_Call SNPs using SAMTOOLS v0.1.19_

```bash
# SNPs and variants are called across all samples simultaneously
# Make list of bam files
ls *.sorted.rmdup.mq20.real.bam > bamfiles.txt

# Call SNPs with SAMTOOLS v.0.1.19
# Filter for SNPs with a base quality score ≥20,
# on reads with a mapping quality ≥20
# BAQ - or realignment process is automatically enabled.
samtools \
   mpileup \
   -ugD \
   -q20 \
   -Q20 \
   -S \
   -d 200 \
   -f CB1.fasta \
   -b bamfiles.txt | \
   bcftools \
      view -vcg -> All-samtools.vcf
   
   # Results: 18.176,225 SNPs, 3,204,056 INDELS, and 21,380,281 total
```

_Call SNPs using ANGSD v0.563_

```bash
# Call SNPs with ANGSD
# Filter for SNPs on reads with a mapping quality ≥20,
# include -baq realignment from SAMTOOLS,
# exclude genotypes with posterior probability <0.99
# LRT = likelihood ratio test cutoff for Hardy-Weinberg equilibrium
angsd \
   -out All-angsd \
   -bam bamfiles.txt \
   -minMapQ 20 \
   -baq \
   -doCounts 1 \
   -GL 1 \
   -uniqueOnly 1 \
   -doGeno 5 \
   -doPost 1 \
   -postCutoff 0.99 \
   -doMajorMinor 4 \
   -doSNP 1 \
   -minLRT 24.0 \
   -doMaf 2 \
   -ref CB1.fasta \
   -nThreads 16
   
   # Results: 19.206,764 SNPs
```

_Get the SNPs overlapping all three algorithms_

```bash
# First convert the ANGSD results to BED format: scaffold\t[pos-1]\tpos
gunzip All-angsd.geno.gz
cut \
   -f1,2 \
   All-angsd.geno | \
   perl -ne '@tmp=split(/\t/,$_); 
      $pos=$tmp[1]-1; 
      print "$tmp[0]\t$pos\t$tmp[1]"' \
      > All-angsd.bed

# Select only SNPs from the SAMTOOLS and GATK results:
vcftools \
   --vcf All-samtools.vcf \
   --remove-indels \
   --recode \
   --recode-INFO-all \
   --out All-samtools.SNPs.vcf
mv All-samtools.SNPs.vcf.recode.vcf All-samtools.SNPs.vcf

vcftools \
   --vcf All-GATK.vcf \
   --remove-indels \
   --recode \
   --recode-INFO-all \
   --out All-GATK.SNPs.vcf
mv All-GATK.SNPs.vcf.recode.vcf All-GATK.SNPs.vcf

# Find pairwise overlapping SNP using BEDTOOLS 2.17.0
intersectBed \
   -a All-GATK.SNPs.vcf \
   -b All-samtools.SNPs.vcf \
   -wa \
   -header > All-GATK-SAMTOOLS.overlap.vcf
intersectBed \
   -a All-GATK.SNPs.vcf \
   -b All-angsd.bed \
   -wa -header > All-GATK-ANGSD.overlap.vcf
intersectBed \
   -a All-samtools.SNPs.vcf \
   -b All-angsd.bed \
   -wa -header > All-SAMTOOLS-ANGSD.overlap.vcf
intersectBed \
   -a All-GATK-SAMTOOLS.overlap.vcf \
   -b All-angsd.bed \
   -wa \
   -header > All-GATK-SAMTOOLS-ANGSD.overlap.vcf
   
   # Results: 15,972,539 overlapping SNPs
```

_Perform stringent SNP filtering to generate a high-quality SNP dataset to use for base quality score recalibration (BQSR)_

```bash
# Download Repetive regions from GENBANK (masking coordinates) for reference
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Camelus_ferus/masking_coordinates.gz
gunzip masking_coordinates.gz
mv masking_coordinates masking_coordinates.bed

# Filter overlapping SNPs from repetitive regions
intersectBed \
   -a All-GATK-SAMTOOLS-ANGSD.overlap.vcf \
   -b masking_coordinates.bed \
   -wa \
   -header > All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.vcf

	# RESULTS
	# All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.vcf: 5,219,307 SNPs
	# A total of 10,753,232 SNPs (67.3%) removed

# Filter SNPs using various GATK Statistics
   # Genotype quality < 10
   # coverage > 750X across all SNPs
   # Quality x depth score < 2
   # Fisher strand bias test > 60
   # Mapping quality < 40
   # Inbreeding coefficient < -0.8 (excess heterozygosity)
   # Mapping quality rank sum test < -12.5
   # Read position rank sum test < -8.0
   # 3 or more SNPs in any 20 bp window
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T VariantFiltration \
   -R CB1.fasta \
   -o All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.vcf \
   --variant All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.vcf \
   --genotypeFilterExpression "GQ < 10.0" \
   --genotypeFilterName "GQ10" \
   --filterExpression "DP > 750 || QD < 2.0 || FS > 60.0 || MQ < 40.0" \
   --filterName "GATKfilter" \
   --filterExpression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
   --filterName "ExcessHet" \
   --filterExpression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
   --filterName "MQRSlt12" \
   --filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
   --filterName "RPRSlt8" \
   --clusterSize 3 \
   --clusterWindowSize 20

# Select variants GATK
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R CB1.fasta \
   -o All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.vcf \
   --variant All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.vcf \
   -ef \
   -env \
   --restrictAllelesTo BIALLELIC 
	# RESULTS
	# All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.vcf:
	# 4,580,891 SNPs

# Find SNPs too close to GATK-called INDELS
# First get indel-only VCF file
vcftools \
   --vcf All-GATK.vcf \
   --keep-only-indels \
   --recode \
   --recode-INFO-all \
   --out All-GATK.indels.vcf
mv All-GATK.indels.vcf.recode.vcf All-GATK.indels.vcf

# Convert to bed file and add 10bp flanking to each side of indels
vcf2bed \
   --deletions \
   --do-not-sort < All-GATK.indels.vcf > All-GATK.dels
vcf2bed \
   --insertions \
   --do-not-sort < All-GATK.indels.vcf > All-GATK.ins
vcf2bed \
   --snvs \
   --do-not-sort < All-GATK.indels.vcf > All-GATK.snvs
bedops \
   --range 10 \
   --everything All-GATK.{dels,ins,snvs} > All-GATK.indels.bed

# Now intersect with SNPs:
intersectBed \
   -v \
   -header \
   -a All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.vcf \
   -b All-GATK.indels.bed \
   > All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.IndelGap10.vcf

	# RESULTS
	# All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.IndelGap10.vcf:
	# 4,423,472 SNPs
```

_Perform BQSR in GATK_

```bash
# BQSR before
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -R CB1.fasta \
   -I Drom439.sorted.rmdup.mq20.real.bam \
   -I Drom795.sorted.rmdup.mq20.real.bam \
   -I Drom796.sorted.rmdup.mq20.real.bam \
   -I Drom797.sorted.rmdup.mq20.real.bam \
   -I Drom800.sorted.rmdup.mq20.real.bam \
   -I Drom802.sorted.rmdup.mq20.real.bam \
   -I Drom806.sorted.rmdup.mq20.real.bam \
   -I Drom816.sorted.rmdup.mq20.real.bam \
   -I Drom820.sorted.rmdup.mq20.real.bam \
   -I DC158.sorted.rmdup.mq20.real.bam \
   -I DC269.sorted.rmdup.mq20.real.bam \
   -I DC399.sorted.rmdup.mq20.real.bam \
   -I DC400.sorted.rmdup.mq20.real.bam \
   -I DC402.sorted.rmdup.mq20.real.bam \
   -I DC408.sorted.rmdup.mq20.real.bam \
   -I DC423.sorted.rmdup.mq20.real.bam \
   -I WC214.sorted.rmdup.mq20.real.bam \
   -I WC216.sorted.rmdup.mq20.real.bam \
   -I WC218.sorted.rmdup.mq20.real.bam \
   -I WC219.sorted.rmdup.mq20.real.bam \
   -I WC220.sorted.rmdup.mq20.real.bam \
   -I WC247.sorted.rmdup.mq20.real.bam \
   -I WC303.sorted.rmdup.mq20.real.bam \
   -I WC304.sorted.rmdup.mq20.real.bam \
   -I WC305.sorted.rmdup.mq20.real.bam \
   -knownSites \
      All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.IndelGap10.vcf \
   --bqsrBAQGapOpenPenalty 30 \
   -o All_BQSR.before.grp \
   -dcov 200 \
   -nct 16

# BQSR After
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -R CB1.fasta \
   -I Drom439.sorted.rmdup.mq20.real.bam \
   -I Drom795.sorted.rmdup.mq20.real.bam \
   -I Drom796.sorted.rmdup.mq20.real.bam \
   -I Drom797.sorted.rmdup.mq20.real.bam \
   -I Drom800.sorted.rmdup.mq20.real.bam \
   -I Drom802.sorted.rmdup.mq20.real.bam \
   -I Drom806.sorted.rmdup.mq20.real.bam \
   -I Drom816.sorted.rmdup.mq20.real.bam \
   -I Drom820.sorted.rmdup.mq20.real.bam \
   -I DC158.sorted.rmdup.mq20.real.bam \
   -I DC269.sorted.rmdup.mq20.real.bam \
   -I DC399.sorted.rmdup.mq20.real.bam \
   -I DC400.sorted.rmdup.mq20.real.bam \
   -I DC402.sorted.rmdup.mq20.real.bam \
   -I DC408.sorted.rmdup.mq20.real.bam \
   -I DC423.sorted.rmdup.mq20.real.bam \
   -I WC214.sorted.rmdup.mq20.real.bam \
   -I WC216.sorted.rmdup.mq20.real.bam \
   -I WC218.sorted.rmdup.mq20.real.bam \
   -I WC219.sorted.rmdup.mq20.real.bam \
   -I WC220.sorted.rmdup.mq20.real.bam \
   -I WC247.sorted.rmdup.mq20.real.bam \
   -I WC303.sorted.rmdup.mq20.real.bam \
   -I WC304.sorted.rmdup.mq20.real.bam \
   -I WC305.sorted.rmdup.mq20.real.bam \
   -knownSites \
      All-GATK-SAMTOOLS-ANGSD.overlap.NoRepeats.filtered.selected.IndelGap10.vcf \
   --bqsrBAQGapOpenPenalty 30 \
   -BQSR All_BQSR.before.grp \
   -o All_BQSR.after.grp \
   -dcov 200 \
   -nct 16

# Make BQSR Plots
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T AnalyzeCovariates \
   -R CB1.fasta \
   -before All_BQSR.before.grp \
   -after All_BQSR.after.grp \
   -plots All-BQSR-plots.pdf

# Print New Recalibrated bam files
while read i
   do
   java -Xmx20g -jar GenomeAnalysisTK.jar \
      -T PrintReads \
      -R CB1.fasta \
      -I ${i}.sorted.rmdup.mq20.real.bam \
      -BQSR All_BQSR.before.grp \
      -o ${i}.sorted.rmdup.mq20.real.recal.bam
   done < IDs.txt
```
__The results of the BQSR can be viewed [Here](./plots/BQSR.pdf).__  

_Perform SNP calling using GATK's HaplotypeCaller on the recalibrated BAM files from above_

```bash
java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R CB1.fasta \
	-I Drom439.sorted.rmdup.mq20.real.recal.bam \
	-I Drom795.sorted.rmdup.mq20.real.recal.bam \
	-I Drom796.sorted.rmdup.mq20.real.recal.bam \
	-I Drom797.sorted.rmdup.mq20.real.recal.bam \
	-I Drom800.sorted.rmdup.mq20.real.recal.bam \
	-I Drom802.sorted.rmdup.mq20.real.recal.bam \
	-I Drom806.sorted.rmdup.mq20.real.recal.bam \
	-I Drom816.sorted.rmdup.mq20.real.recal.bam \
	-I Drom820.sorted.rmdup.mq20.real.recal.bam \
	-I DC158.sorted.rmdup.mq20.real.recal.bam \
	-I DC269.sorted.rmdup.mq20.real.recal.bam \
	-I DC399.sorted.rmdup.mq20.real.recal.bam \
	-I DC400.sorted.rmdup.mq20.real.recal.bam \
	-I DC402.sorted.rmdup.mq20.real.recal.bam \
	-I DC408.sorted.rmdup.mq20.real.recal.bam \
	-I DC423.sorted.rmdup.mq20.real.recal.bam \
	-I WC214.sorted.rmdup.mq20.real.recal.bam \
	-I WC216.sorted.rmdup.mq20.real.recal.bam \
	-I WC218.sorted.rmdup.mq20.real.recal.bam \
	-I WC219.sorted.rmdup.mq20.real.recal.bam \
	-I WC220.sorted.rmdup.mq20.real.recal.bam \
	-I WC247.sorted.rmdup.mq20.real.recal.bam \
	-I WC303.sorted.rmdup.mq20.real.recal.bam \
	-I WC304.sorted.rmdup.mq20.real.recal.bam \
	-I WC305.sorted.rmdup.mq20.real.recal.bam \
	-nct 4 \
	--genotyping_mode DISCOVERY \
	-stand_emit_conf 20 \
	-stand_call_conf 20 \
	-dcov 200 \
	-o All-GATKHC-2.vcf
	
	# Result: 19,360,771 SNPs
```

_Get the depth of coverage from the BAM files in order to use the coverage levels for filtering SNPs_

```bash
# Make a list of the bam files
ls *.sorted.rmdup.mq20.real.recal.bam > bamfiles.txt

# Get coverage using either Samtools or Qualimap
# Then analyze as a separate file using stats software (e.g., R)
samtools \
   depth \
   -f bamfiles.txt > coverage.txt

# Qualimap v2.0.1
while read bam
do
qualimap bamqc \
	-bam $bam \
	-nr 2000 \
	-nt 16 \
	-nw 400 \
	-oc $bam.coverage \
	-outdir QUALIMAP/$bam \
	-outfile $bam.pdf \
	-outformat pdf \
	--java-mem-size=3500M
done < bamfiles.txt
   
	# Results:
	# Drom439 = 14.68X +- 16.55X
	# Drom795 = 11.65X +- 13.27X
	# Drom796 = 14.6X +- 15.91X
	# Drom797 = 14.11X +- 13.88X
	# Drom800 = 14.44X +- 16.65X
	# Drom802 = 14.9X +- 15.5X
	# Drom806 = 9.69X +- 11.6X
	# Drom816 = 10.46X +- 11.67X
	# Drom820 = 9.9X +- 10.76X
	# DC158 = 13.78X +- 16.53X
	# DC269 = 14.34X +- 15.78X
	# DC399 = 13.82X +- 16.17X
	# DC400 = 14.61X +- 16.61X
	# DC402 = 14.91X +- 17.09X
	# DC408 = 15.18X +- 17.71X
	# DC423 = 14.58X +- 17.01X
	# WC214 = 14.89X +- 17.71X
	# WC216 = 13.51X +- 16.78X
	# WC218 = 14.69X +- 16.58X
	# WC219 = 14.53X +- 17.92X
	# WC220 = 15.4X +- 16.87X
	# WC247 = 14.41X +- 16.03X
	# WC303 = 14.57X +- 18.89X
	# WC304 = 15.14X +- 17.52X
	# WC305 = 14.37X +- 17.63X

# mean of mean depths = 13.94X ; mean of SD = 16.17X
# Set Max Depth to (mean+1SD) ~ 30X per individual
```

_Filter the SNPs identified above using GATK v3.1-1_
```bash
# First separate SNPs, INDELS/MIXED for separate processing (MIXED= combination of SNP and INDEL at same site)
java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R CB1.fasta \
	--variant All-GATKHC-2.vcf \
	-selectType SNP \
	-o All.raw.SNPs.vcf
	# Results:  16,557,930 SNPs

java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R CB1.fasta \
	--variant All-GATKHC-2.vcf \
	-selectType INDEL \
	-o All.raw.INDELs.vcf
	# Results: 2,721,748 INDELs

java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R CB1.fasta \
	--variant All-GATKHC-2.vcf \
	-selectType MIXED \
	-selectType MNP \
	-o All.raw.MIXED.vcf
	# Results: 81,093 MIXED

# Filter SNP File using GATK statistics
   # Genotype quality < 13
   # coverage > 933X across all SNPs
   # Quality x depth score < 2
   # Fisher strand bias test > 60
   # Mapping quality < 40
   # Inbreeding coefficient < -0.8 (excess heterozygosity)
   # Mapping quality rank sum test < -12.5
   # Read position rank sum test < -8.0
   # 3 or more SNPs in any 20 bp window

# SNP Cluster: if three SNPs are in a 20 bp window.
java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R CB1.fasta \
	-o All.filtered.SNPs.vcf \
	--variant All.raw.SNPs.vcf \
	--filterExpression "DP > 933 || DP < 5" \
	--filterName "DP" \
	--filterExpression "FS > 60.0" \
	--filterName "FS" \
	--filterExpression "MQ < 40.0" \
	--filterName "MQ" \
	--filterExpression "vc.hasAttribute('QD') && QD < 2.0" \
	--filterName "QD" \
	--filterExpression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
	--filterName "F" \
	--filterExpression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
	--filterName "MQRS" \
	--filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
	--filterName "RPRS" \
	--genotypeFilterExpression "GQ < 13.0" \
	--genotypeFilterName "GQ13" \
	--clusterSize 3 \
	--clusterWindowSize 20

# Filter Indels
java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R CB1.fasta \
	-o All.filtered.INDELs.vcf \
	--variant All.raw.INDELs.vcf \
	--filterExpression "DP > 933 || DP < 5" \
	--filterName "DP" \
	--filterExpression "vc.hasAttribute('QD') && QD < 2.0" \
	--filterName "QD" \
	--filterExpression "FS > 200.0" \
	--filterName "FS" \
	--filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
	--filterName "RPRS"

# Filter MIXED variants
java -Xmx20g -jar GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R CB1.fasta \
	-o All.filtered.MIXED.vcf \
	--variant All.raw.MIXED.vcf \
	--filterExpression "DP > 933 || DP < 5" \
	--filterName "DP" \
	--filterExpression "vc.hasAttribute('QD') && QD < 2.0" \
	--filterName "QD" \
	--filterExpression "FS > 200.0" \
	--filterName "FS" \
	--filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
	--filterName "RPRS"

```

# Variant Quality Score Recalibration (VQSR) using GATK

_Make a bed file of the XY scaffolds to exclude and the masking coordinates (repetitive elements)
```bash
# The file "masking_coordinates" was downloaded earlier from the CB1 reference genome (see code above)
cat \
	<(echo "# BED header") \
	<(cat masking_coordinates | cut -f1-3) > exclude.bed.tmp

# Make all genomic windows based on the reference genome index (.fai) created earlier using Samtools
bedtools makewindows \
   -n 1 \
   -g CB1.fasta.fai > genome.bed

# Find matching XY scaffolds from the genome.bed file and add to the list of excluded regions
grep -f XY.exclude genome.bed >> exclude.bed.tmp

# Clean up bed file (sort and merge regions
cat \
	<(echo "# BED header") \
	<(bedtools sort -header -i exclude.bed.tmp | bedtools merge -i /dev/stdin) > exclude.bed
```

_Additional filtering of variants usinf VCFTOOLS_
```bash
# Finish cleaning of SNPs using VCFTOOLS v0.1.12b
   # Remove all filtered variants and genotypes from GATK
   # minimum allele count must be ≥2
   # Remove variants from regions in bed file
   # Remove variants with a HWE p-value < 0.001 (not in HWE)
   # Genotype depth (individual) must be ≥5 and ≤30
vcftools --vcf All.filtered.SNPs.vcf \
	--remove-filtered-all \
	--remove-filtered-geno-all \
	--exclude-bed exclude.bed \
	--mac 2 \
	--max-missing-count 5 \
	--hwe 0.0001 \
	--minDP 5 \
	--maxDP 30 \
	--recode \
	--recode-INFO-all \
	--out VQSR.training.vcf
	# Results:
	# Outputting VCF file...
	# 	Read 2750974 BED file entries.
	# After filtering, kept 4869698 out of a possible 16557930 Sites
```
