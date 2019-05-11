# Homogeneity Tests for Selection in Camels
The homogeneity test is a test of excess polymorphism vs divergence and first described by [Liu et al 2014](http://dx.doi.org/10.1016/j.cell.2014.03.054) when resequencing polar bear genomes. The test is performed on a gene-by-gene basis. The homogeneity test can be 'polarized' to determine which taxa most likely experienced the loss of polymorphism using a genome-wide HKA test.  Please see our manuscript or the reference above for further details on the methodology.  The procedure can be boiled down into three steps:
1. Estimate the ancestral allele at each site using the alpaca genome
2. Perform the homogeneity test to examine the intraspecific (polymorphism) and interspecific (divergence) genetic diversity
3. Perform an HKA test for each lineage of the test.

## Step 1: Make a list of SNPs polymorphic within each camel species_
The files [Drom.txt](./Data/Drom.txt), [DC.txt](./Data/DC.txt), and [WC.txt](./Data/WC.txt) are simply lists of the individual IDs.  The SNPs were annotated using SNPEFF v.

```bash
# Dromedaries (Drom)
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep Drom.txt \
   --maf 0.01 \
   --max-missing-count 2 \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s Drom.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > Drom.EFF.vcf 2> Drom.out
   # Results: After filtering, kept 2656809 out of a possible 10819573 Sites
   
# Wild camels (WC)
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep WC.txt \
   --maf 0.01 \
   --max-missing-count 2 \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s WC.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > WC.EFF.vcf 2> WC.out
   # Results: After filtering, kept 3899346 out of a possible 10819573 Sites

# Domestic Bactrian camels (DC)
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep DC.txt \
   --maf 0.01 \
   --max-missing-count 1 \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s DC.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > DC.EFF.vcf 2> DC.out
   # After filtering, kept 5000559 out of a possible 10819573 Sites
```

## Step 2: Make a list of SNPs fixed (F<sub>ST</sub> = 1) between camel species
_Drom vs WC groups_
```bash
# DROM vs WC Fst calculations
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --weir-fst-pop Drom.txt \
   --weir-fst-pop WC.txt \
   --stdout | \
   sed '1d' | \
   perl -ne 'chomp;@a=split(/\t/,$_);if($a[2]==1){print "$a[0]\t$a[1]\n";}' \
   > Drom-WC.Fst.pos
   # Results: 3410975 SNPs
```

_DC vs WC groups_
```bash
# DC vs WC Fst calculations
vcftools 
   --vcf All.SNPs.filtered.vcf \
   --weir-fst-pop DC.txt \
   --weir-fst-pop WC.txt \
   --stdout | \
   sed '1d' | \
   perl -ne 'chomp;@a=split(/\t/,$_);if($a[2]==1){print "$a[0]\t$a[1]\n";}' \
   > DC-WC.Fst.pos
	# Results: 18801 SNPs
```

## Step 3: Estimate ancestral alleles at each site using the alpaca genome as an outgroup

_Map PE, Illumina alpaca genome reads from [Wu et al. 2014](https://doi.org/10.1038/ncomms6188) to the camel reference_
The SRA accession numbers of the reads were:
1. [SRR1552599](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552599)
2. [SRR1552605](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552605)
3. [SRR1552606](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552606)
4. [SRR1552607](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552607)
5. [SRR1552610](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552610)

```bash
# Download and map reads for each accession above (only SRR1552599 is shown as an example)
lib="SRR1552599"

# Download reads in fastq format using SRATOOLKIT
fastq-dump \
   -v -v \
   --split-files \
   $lib

# Trim reads using POPOOLATION v1.2.2 as for the camel reads
perl trim-fastq.pl \
   --input1 ${lib}_1.fastq \
   --input2 ${lib}_2.fastq \
   --output ${lib}.trimmed \
   --quality-threshold 20 \
   -fastq-type sanger \
   --min-length 50 \
   > ${lib}.stats

# Set some values for the 'Read Groups tag in samtools
rg="@RG\tID:CB1\tPL:illumina\tPU:1\tLB:${lib}\tSM:${lib}\tCN:BGI"

# Map reads with BWA v0.6.2
bwa sampe \
   -r $rg \
   CB1.fasta \
   <(bwa aln -n 0.01 -o 1 -e 12 -d 12 -l 32 -t 16 CB1.fasta ${lib}.trimmed_1.fq) \
   <(bwa aln -n 0.01 -o 1 -e 12 -d 12 -l 32 -t 16 CB1.fasta ${lib}.trimmed_2.fq) \
   ${lib}.trimmed_1.fq \
   ${lib}.trimmed_2.fq | \
   samtools view -Shu - | \
   samtools sort - ${lib}.sorted

# Remove duplicate reads, filter for MQ > 20, properly paired and mapped reads, convert to BAM
samtools \
   rmdup \
   -s ${lib}.sorted.bam - | \
   samtools \
      view \
      -bh \
      -q 20 \
      -f 0x0002 \
      -F 0x0004 \
      -F 0x0008 - > ${lib}.sorted.rmdup.mq20.bam

# Index bam file
samtools index ${lib}.sorted.rmdup.mq20.bam

# Get mapping statistics
samtools \
   stats ${lib}.sorted.rmdup.mq20.bam > ${lib}.bamstats
```

_Merge the BAM files for each SRR file_
```bash
# Make a header from the last bam file
samtools \
   view \
   -H ${lib}.sorted.rmdup.mq20.bam > header.txt

# Make a list of the 5 BAM files
ls *.sorted.rmdup.mq20.bam > bams.list

# Merge the bam files using SAMTOOLS v1.3
samtools \
   merge \
   -h header.txt \
   -b bams.list \
   Vpacos.merged.bam

# Index the alpaca merged bam file
samtools index Vpacos.merged.bam

# get overall mapping statistics
samtools stats \
   Vpacos.merged.bam > Vpacos.merged.bamstats
```

_Use ANGSD v0.563 to make allele counts at each Fixed SNP site_
```bash
# Index the sites for [Drom vs WC] and [DC vs WC]
angsd sites index Drom-WC.Fst.pos
angsd sites index DC-WC.Fst.pos

# Run ANGSD to make base counts at each site (SNP)
# Drom vs WC
angsd \
   -out Vpacos.counts.Drom-WC \
   -i Vpacos.merged.bam \
   -doCounts 1 \
   -dumpCounts 3 \
   -baq 1 \
   -nThreads 16 \
   -ref CB1.fasta \
   -sites ../Drom-WC.Fst.pos

# DC vs WC
angsd \
   -out Vpacos.counts \
   -i Vpacos.merged.bam \
   -doCounts 1 \
   -dumpCounts 3 \
   -baq 1 \
   -nThreads 16 \
   -ref CB1.fasta \
   -sites DC-WC.Fst.pos
```





# Subset the VCF file for these SNPs and annotate them using SNPEFF
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --positions Drom-WC.Fst.pos \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s DROM-WC.fixed.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > DROM-WC.fixed.EFF.vcf 2> DROM-WC.fixed.out
   # Results: After filtering, kept 3410975 out of a possible 10819573 Sites
```











   # Run code at bottom to get new .pos file
vcftools --vcf /wrk/rfitak/SNP-ANALYSIS/VQSR/All.SNPs.filtered.vcf \
	--positions ALPACA/DC-WC-Vpacos.fixed.pos \
	--recode \
	--stdout | \
	java -Xmx4g -jar /wrk/rfitak/SOFTWARE/snpEff/snpEff.jar ann \
	-v \
	-i vcf \
	-o vcf \
	-s DC-WC-Vpacos.fixed.stats.html \
	-canon \
	-onlyProtein \
	camelus_ferus \
	- > DC-WC-Vpacos.fixed.EFF.vcf 2> DC-WC-Vpacos.fixed.out
	# Results: After filtering, kept 6745 out of a possible 10819573 Sites
vcftools --vcf /wrk/rfitak/SNP-ANALYSIS/VQSR/All.SNPs.filtered.vcf \
	--positions ALPACA/WC-DC-Vpacos.fixed.pos \
	--recode \
	--stdout | \
	java -Xmx4g -jar /wrk/rfitak/SOFTWARE/snpEff/snpEff.jar ann \
	-v \
	-i vcf \
	-o vcf \
	-s WC-DC-Vpacos.fixed.stats.html \
	-canon \
	-onlyProtein \
	camelus_ferus \
	- > WC-DC-Vpacos.fixed.EFF.vcf 2> WC-DC-Vpacos.fixed.out
	# Results: After filtering, kept 11567 out of a possible 10819573 Sites
vcftools --vcf /wrk/rfitak/SNP-ANALYSIS/VQSR/All.SNPs.filtered.vcf \
	--positions ALPACA/Drom-WC-Vpacos.fixed.pos \
	--recode \
	--stdout | \
	java -Xmx4g -jar /wrk/rfitak/SOFTWARE/snpEff/snpEff.jar ann \
	-v \
	-i vcf \
	-o vcf \
	-s Drom-WC-Vpacos.fixed.stats.html \
	-canon \
	-onlyProtein \
	camelus_ferus \
	- > Drom-WC-Vpacos.fixed.EFF.vcf 2> Drom-WC-Vpacos.fixed.out
	# Results: After filtering, kept 1837134 out of a possible 10819573 Sites
vcftools --vcf /wrk/rfitak/SNP-ANALYSIS/VQSR/All.SNPs.filtered.vcf \
	--positions ALPACA/WC-Drom-Vpacos.fixed.pos \
	--recode \
	--stdout | \
	java -Xmx4g -jar /wrk/rfitak/SOFTWARE/snpEff/snpEff.jar ann \
	-v \
	-i vcf \
	-o vcf \
	-s WC-Drom-Vpacos.fixed.stats.html \
	-canon \
	-onlyProtein \
	camelus_ferus \
	- > WC-Drom-Vpacos.fixed.EFF.vcf 2> WC-Drom-Vpacos.fixed.out
	# Results: After filtering, kept 1523950 out of a possible 10819573 Sites
```
