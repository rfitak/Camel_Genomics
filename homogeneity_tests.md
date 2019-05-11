# Homogeneity Tests for Selection in Camels
The homogeneity test is a test of excess polymorphism vs divergence and first described by [Liu et al 2014](http://dx.doi.org/10.1016/j.cell.2014.03.054) when resequencing polar bear genomes. The test is performed on a gene-by-gene basis. The homogeneity test can be 'polarized' to determine which taxa most likely experienced the loss of polymorphism using a genome-wide HKA test.  Please see our manuscript or the reference above for further details on the methodology.  The procedure can be boiled down into three steps:
1. Estimate the ancestral allele at each site using the alpaca genome
2. Perform the homogeneity test to examine the intraspecific (polymorphism) and interspecific (divergence) genetic diversity
3. Perform an HKA test for each lineage of the test.

## Step 1: Make a list of SNPs polymorphic within each camel species
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
   -out Vpacos.counts.DC-WC \
   -i Vpacos.merged.bam \
   -doCounts 1 \
   -dumpCounts 3 \
   -baq 1 \
   -nThreads 16 \
   -ref CB1.fasta \
   -sites DC-WC.Fst.pos
```

_From the counts files, the following R code selects the ancestral allele_
```R
# DROM vs WC
# Read in the ANGSD output files
a = read.table(gzfile("Vpacos.counts.Drom-WC.pos.gz"), header = T)
b = read.table(gzfile("Vpacos.counts.Drom-WC.counts.gz"), header = T)

# Combine output into a matrix
df = cbind(a, b)

# Retain sites with a depth of at least 2
df2 = subset(df, totDepth >= 2)
   # Started with 3410975 positions
      # 3332844 with at least 1x coverage
      # 3301577 with >=2x coverage

# Make a list of bases to count (no 'N')
bases = c("A", "C", "G", "T")

# For each SNP, select the base with the most counts
# If there are ties, randomly select one of the tied bases
Vpacos = vector()
for (i in 1:nrow(df2)){
   mx = which(df2[i, 4:7] == max(df2[i, 4:7]))
   if (length(mx) == 1){
      anc = mx
   } else {
      anc = sample(mx, 1)
   }
   Vpacos = c(Vpacos, bases[anc])
   print(i)
}

# Append a column of the ancestral base selected
df2 = cbind(df2, Vpacos)

# Write the output to a new table
write.table(df2, file = "Vpacos.Drom-WC.anc.tsv", sep = "\t", row.names = F, quote = F)

# DC vs WC
# Read in the ANGSD output files
a = read.table(gzfile("Vpacos.counts.DC-WC.pos.gz"), header = T)
b = read.table(gzfile("Vpacos.counts.DC-WC.counts.gz"), header = T)

# Combine output into a matrix
df = cbind(a, b)

# Retain sites with a depth of at least 2
df2 = subset(df, totDepth >= 2)
   # Started with 18801 positions
      # 18251 with at least 1x coverage
      # 18001 with >=2x coverage

# Make a list of bases to count (no 'N')
bases = c("A", "C", "G", "T")

# For each SNP, select the base with the most counts
# If there are ties, randomly select one of the tied bases
Vpacos = vector()
for (i in 1:nrow(df2)){
   mx = which(df2[i, 4:7] == max(df2[i, 4:7]))
   if (length(mx) == 1){
      anc = mx
   } else {
      anc = sample(mx, 1)
   }
   Vpacos = c(Vpacos, bases[anc])
}

# Append a column of the ancestral base selected
df2 = cbind(df2, Vpacos)

# Write the output to a new table
write.table(df2, file = "Vpacos.DC-WC.anc.tsv", sep = "\t", row.names = F, quote = F)
```

_Get the alleles from WC and Drom to append to the table_  
see the [anc.pl](./Data/anc.pl) script.
```bash
# WC - append column to table titled "WC_Allele"
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep WC.txt \
   --positions Vpacos.Drom-WC.anc.tsv \
   --stdout \
   --counts | \
   sed '1d' | \
   ./anc.pl | \
   paste Vpacos.Drom-WC.anc.tsv <(cat <(echo "WC_Allele") -) > tmp
mv tmp Vpacos-WC-Drom.pos

# DROM - append column to table titled "DROM_Allele"
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep Drom.txt \
   --positions Vpacos.Drom-WC.anc.tsv \
   --stdout \
   --counts | \
   sed '1d' | \
   ./anc.pl | \
   paste Vpacos-WC-Drom.pos <(cat <(echo "DROM_Allele") -) > tmp
mv tmp Vpacos-WC-Drom.pos
```

_Make a list of sites fixed between either DROM vs \[WC+alpaca\] or  WC vs \[DROM+alpaca\]_
```bash
# Sites DROM vs [WC+alpaca]
sed '1d' Vpacos-WC-Drom.pos | \
   perl -ne 'chomp; @a=split(/\t/,$_);if($a[9] ne $a[8] && $a[9] ne $a[7]){print "$_\n"}' | \
   cut -f1,2 > Drom-WC-Vpacos.fixed.pos
      # Results
      # 1837134 SNPs to process
      
# Sites WC vs [DROM+alpaca]
sed '1d' Vpacos-WC-Drom.pos | \
   perl -ne 'chomp; @a=split(/\t/,$_);if($a[8] ne $a[7] && $a[8] ne $a[9]){print "$_\n"}' | \
   cut -f1,2 > WC-Drom-Vpacos.fixed.pos
      # Results
      # 1523950 SNPs to process
```

_Get the alleles from WC and DC to append to the table_
```bash
# WC - append column to table titled "WC_Allele"
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep WC.txt \
   --positions Vpacos.DC-WC.anc.tsv \
   --stdout \
   --counts | \
   sed '1d' | \
   ./anc.pl | \
   paste Vpacos.DC-WC.anc.tsv <(cat <(echo "WC_Allele") -) > tmp
mv tmp Vpacos-WC-DC.pos

# DC - append column to table titled "DC_Allele"
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --keep DC.txt \
   --positions Vpacos.DC-WC.anc.tsv \
   --stdout \
   --counts | \
   sed '1d' | \
   ./anc.pl | \
   paste Vpacos-WC-DC.pos <(cat <(echo "DC_Allele") -) > tmp
mv tmp Vpacos-WC-DC.pos
```

_Make a list of sites fixed between either DC vs \[WC+alpaca\] or  WC vs \[DC+alpaca\]_
```bash
# Sites DC vs [WC+alpaca]
sed '1d' Vpacos-WC-DC.pos | \
   perl -ne 'chomp; @a=split(/\t/,$_);if($a[9] ne $a[8] && $a[9] ne $a[7]){print "$_\n"}' | \
   cut -f1,2 > DC-WC-Vpacos.fixed.pos
      # Results
      # 6745 SNPs to process
      
# Sites WC vs [DC+alpaca]
sed '1d' Vpacos-WC-DC.pos | \
   perl -ne 'chomp; @a=split(/\t/,$_);if($a[8] ne $a[7] && $a[8] ne $a[9]){print "$_\n"}' | \
   cut -f1,2 > WC-DC-Vpacos.fixed.pos
      # Results
      # 11567 SNPs to process
```

_Make subsets of the VCF file for these Fixed SNPs and annotate using SNPEFF_
```bash
# Drom vs WC
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --positions Drom-WC-Vpacos.fixed.pos \
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

# WC vs Drom
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --positions WC-Drom-Vpacos.fixed.pos \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s WC-Drom-Vpacos.fixed.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > WC-Drom-Vpacos.fixed.EFF.vcf 2> WC-Drom-Vpacos.fixed.out
   # Results: After filtering, kept 1523950 out of a possible 10819573 Sites

# DC vs WC
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --positions DC-WC-Vpacos.fixed.pos \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s DC-WC-Vpacos.fixed.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > DC-WC-Vpacos.fixed.EFF.vcf 2> DC-WC-Vpacos.fixed.out
   # Results: After filtering, kept 6745 out of a possible 10819573 Sites

# WC vs DC
vcftools \
   --vcf All.SNPs.filtered.vcf \
   --positions ALPACA/WC-DC-Vpacos.fixed.pos \
   --recode \
   --stdout | \
   java -Xmx4g -jar snpEff.jar ann \
   -v \
   -i vcf \
   -o vcf \
   -s WC-DC-Vpacos.fixed.stats.html \
   -canon \
   -onlyProtein \
   camelus_ferus \
   - > WC-DC-Vpacos.fixed.EFF.vcf 2> WC-DC-Vpacos.fixed.out
   # Results: After filtering, kept 11567 out of a possible 10819573 Sites   
```

## Step 4: Perform homogeneity and HKA tests
This following code basically builds a big table of all the results and statistical tests which can be later filtered for the desired criteria of interest (e.g., using EXCEL)

_Get a list of all transcripts annotated in SNPEFF_
```bash
sed '1d' All.SNPs.filtered.stats.genes.txt | \
   sed '1d' | \
   cut -f3 > transcripts.list
   # Results: 17,912 transcripts
```

_Build some custom R functions to use in the analysis_
```R
# Build custom functions
fisher <- function(x){
	x=matrix(as.integer(x[2:5]),nrow=2, byrow=T)
	z=fisher.test(x)
	return(z$p.value)
}
HKA_AC <- function(x){
	x=c(x[2], A.hka, x[4], C.hka)
	x=matrix(as.integer(x),nrow=2, byrow=T)
	z=fisher.test(x)
	return(z$p.value)
}
HKA_BD <- function(x){
	x=c(x[3], B.hka, x[5], D.hka)
	x=matrix(as.integer(x),nrow=2, byrow=T)
	z=fisher.test(x)
	return(z$p.value)
}
ratio <- function(x){
	x=as.integer(x[2:5])
	z=(x[3]/x[1]) / (x[4]/x[2])
	return(z)
}
```

_For DROM vs WC, build a table of SNP sites to use_  
In each case, the "bases_affected_EXON" column from SNPEFF is counted.
```bash
# Drom vs WC
c=1
while read line
   do
   # A = count number of polymorphic sites in the dromedary samples
   A=$(grep "$line" Drom.stats.genes.txt | cut -f35)
   if [ "$A" = "" ]
      then
      A="0"
   fi
   # B = count number of polymorphic sites in the wild camel samples
   B=$(grep "$line" WC.stats.genes.txt | cut -f36)
   if [ "$B" = "" ]
      then
      B="0"
   fi
   # C = count number of fixed differences between dromedaries and both wild camels and the alpaca genome sequence.
   C=$(grep "$line" Drom-WC-Vpacos.fixed.stats.genes.txt | cut -f32)
   if [ "$C" = "" ]
      then
      C="0"
   fi
   # D = count number of fixed differences between wild camels and both dromedaries and the alpaca genome sequence.
   D=$(grep "$line" WC-Drom-Vpacos.fixed.stats.genes.txt | cut -f33)
   if [ "$D" = "" ]
      then
      D="0"
   fi
   echo -e "$line\t$A\t$B\t$C\t$D" >> Drom-WC.tbl
   echo "Finished transcript $c"
   c=$(( $c + 1 ))
done < transcripts.list
```

_For DC vs WC, build a table of SNP sites to use_  
In each case, the "bases_affected_EXON" column from SNPEFF is counted.

```bash
c=1
while read line
do
   # A = count number of polymorphic sites in the DC samples
   A=$(grep "$line" DC.stats.genes.txt | cut -f35)
   if [ "$A" = "" ]
      then
      A="0"
   fi
   # B = count number of polymorphic sites in the wild camel samples
   B=$(grep "$line" WC.stats.genes.txt | cut -f36)
   if [ "$B" = "" ]
      then
      B="0"
   fi
   # C = count number of fixed differences between DC and both wild camels and the alpaca genome sequence.
   C=$(grep "$line" DC-WC-Vpacos.fixed.stats.genes.txt | cut -f22)
   if [ "$C" = "" ]
      then
      C="0"
   fi
   # D = count number of fixed differences between wild camels and both DC and the alpaca genome sequence.
   D=$(grep "$line" WC-DC-Vpacos.fixed.stats.genes.txt | cut -f25)
   if [ "$D" = "" ]
      then
      D="0"
   fi
   echo -e "$line\t$A\t$B\t$C\t$D" >> DC-WC.tbl
   echo "Finished transcript $c"
   c=$(( $c + 1 ))
done < transcripts.list
```
