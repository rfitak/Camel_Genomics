# Homogeneity Tests for Selection in Camels
The homogeneity test is a test of excess polymorphism vs divergence and first described by [Liu et al 2014](http://dx.doi.org/10.1016/j.cell.2014.03.054) when resequencing polar bear genomes. The test is performed on a gene-by-gene basis. The homogeneity test can be 'polarized' to determine which taxa most likely experienced the loss of polymorphism using a genome-wide HKA test.  Please see our manuscript or the reference above for further details on the methodology.  The procedure can be boiled down into three steps:
1. Estimate the ancestral allele at each site using the alpaca genome
2. Perform the homogeneity test to examine the intraspecific (polymorphism) and interspecific (divergence) genetic diversity
3. Perform an HKA test for each lineage of the test.

## Step 1: Estimate ancestral alleles at each site using the alpaca genome as an outgroup

_Map PE, Illumina alpaca genome reads from [Wu et al. 2014](https://doi.org/10.1038/ncomms6188) to the camel reference_
The SRA accession numbers of the reads were:
1. [SRR1552599](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552599)
2. [SRR1552605](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552605)
3. [SRR1552606](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552606)
4. [SRR1552607](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552607)
5. [SRR1552610](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1552610)

```bash


cd /wrk/rfitak/HOMOGENEITY_TESTS/ALPACA
HOME="/wrk/rfitak"

lib=$(sed -n "$SLURM_ARRAY_TASK_ID"p SRR.list)

./fastq-dump \
   -v -v \
   --split-files \
   $lib

perl /wrk/rfitak/SOFTWARE/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
   --input1 ${lib}_1.fastq \
   --input2 ${lib}_2.fastq \
   --output ${lib}.trimmed \
   --quality-threshold 20 \
   -fastq-type sanger \
   --min-length 50 \
   > ${lib}.stats

bwa=/wrk/rfitak/DROM_ALL/DROM-MAPPING/bwa-0.6.2/bwa
samtools=/wrk/rfitak/SOFTWARE/EXECUTABLES/samtools
rg="@RG\tID:CB1\tPL:illumina\tPU:1\tLB:${lib}\tSM:${lib}\tCN:BGI"

$bwa sampe \
-r $rg \
CB1.fasta \
<($bwa aln -n 0.01 -o 1 -e 12 -d 12 -l 32 -t 16 CB1.fasta ${lib}.trimmed_1.fq) \
<($bwa aln -n 0.01 -o 1 -e 12 -d 12 -l 32 -t 16 CB1.fasta ${lib}.trimmed_2.fq) \
${lib}.trimmed_1.fq \
${lib}.trimmed_2.fq | \
$samtools view -Shu - | \
$samtools sort - ${lib}.sorted

$samtools rmdup -s ${lib}.sorted.bam - | \
$samtools view -bh -q 20 -f 0x0002 -F 0x0004 -F 0x0008 - > ${lib}.sorted.rmdup.mq20.bam
$samtools index ${lib}.sorted.rmdup.mq20.bam
/homeappl/appl_taito/bio/samtools/1.3/bin/samtools stats ${lib}.sorted.rmdup.mq20.bam > ${lib}.bamstats
$samtools view -H ${lib}.sorted.rmdup.mq20.bam > header.txt
ls *.sorted.rmdup.mq20.bam > bams.list
$samtools13 merge \
   -h header.txt \
   -b bams.list \
   Vpacos.merged.bam

$samtools13 index Vpacos.merged.bam

$samtools13 stats \
   Vpacos.merged.bam > Vpacos.merged.bamstats

# Call Ancestral alleles
angsd/angsd \
   sites index ../DC-WC.Fst.pos
angsd/angsd \
   sites index ../Drom-WC.Fst.pos

angsd/angsd \
   -out Vpacos.counts \
   -i Vpacos.merged.bam \
   -doCounts 1 \
   -dumpCounts 3 \
   -baq 1 \
   -nThreads 16 \
   -ref CB1.fasta \
   -sites ../DC-WC.Fst.pos
angsd/angsd \
   -out Vpacos.counts.Drom-WC \
   -i Vpacos.merged.bam \
   -doCounts 1 \
   -dumpCounts 3 \
   -baq 1 \
   -nThreads 16 \
   -ref CB1.fasta \
   -sites ../Drom-WC.Fst.pos

```
