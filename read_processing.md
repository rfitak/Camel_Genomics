# Read Processing and Alignment
_Read trimming using POPOOLATION 1.2.2_

```bash
# Make a list of the sample IDs
echo "DC158 DC269 DC399 DC400 DC402 DC408 DC423 \
   Drom802 Drom439 Drom795 Drom796 Drom797 Drom800 Drom806 Drom816 Drom820 \
   WC214 WC216 WC218 WC219 WC220 WC247 WC303 WC304 WC305" | \
   tr -s " " | tr " " "\n" > IDs.txt
   
# Trim all sequences using POPOOLATION v1.2.2 perl script trim-fastq.pl
while read i
   do
   perl trim-fastq.pl \
      --input1 ${i}_1.fq \
      --input1 ${i}_2.fq \
      --output ${i}trim \
      --quality-threshold 20 \
      --min-length 50 > ${i}_sumstat.txt
   done < IDs.txt
```

_Download and prepare the reference genome_

```bash
# Download and unpack reference genome (Genbank Accession )
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Camelus_ferus/CHR_Un/cfe_ref_CB1_chrUn.fa.gz
gunzip cfe_ref_CB1_chrUn.fa.gz
mv cfe_ref_CB1_chrUn.fa CB1.fasta

# Index the reference genome using BWA 0.6.2
bwa \
   index \
   -a bwtsw \
   CB1.fasta

# Create sequence dictionary using PICARD v1.89
java -jar CreateSequenceDictionary.jar \
   R=CB1.fasta \
   O=CB1.dict

# Index the reference FASTA using SAMTOOLS v0.1.19
samtools faidx CB1.fasta
```
_Map the trimmed reads to the reference_

```bash
# Only individual DC158 is shown as an example
name=DC158

# Set $trim1 and $trim2 to the trimmed reads for that individual
trim1=${name}trim_1
trim2=${name}trim_2

# Set the read group identified for the output BAM file

rg='@RG\tID:1\tPL:illumina\tPU:025752-90\tLB:DC158\tSM:DC158\tCN:BGI'

# Map reads then sort bam file
bwa sampe \
   -r $rg \
   CB1.fasta \
   <(bwa aln -n 0.01 -o 1 -e 12 -d 12 -l 32 -t 16 -I CB1.fasta $trim1 &) \
   <(bwa aln -n 0.01 -o 1 -e 12 -d 12 -l 32 -t 16 -I CB1.fasta $trim2 &) \
   $trim1 \
   $trim2 |\
   samtools view -Shu - | \
   samtools sort - $name.sorted
```

_Clean the BAM files (reads aligned to reference)_

```bash
# Use PICARD to flag and or remove duplicate reads
java -Xmx20g -Dsnappy.disable=true -jar MarkDuplicates.jar \
   REMOVE_DUPLICATES=true \
   I=$name.sorted.bam \
   O=$name.sorted.rmdup.bam \
   M=$name.picmetrics.txt \
   VALIDATION_STRINGENCY=SILENT \
   MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

# Use SAMTOOLS to also remove duplicate reads,
# filter the resulting BAM file for reads with mapping quality of 20,
# and properly paired and mapped reads
samtools \
   rmdup \
   -s $name.sorted.rmdup.bam - | \
   samtools \
      view \
      -bh \
      -q 20 \
      -f 0x0002 \
      -F 0x0004 \
      -F 0x0008 - > $name.sorted.rmdup.mq20.bam

# Index the final bam file
samtools index $name.sorted.rmdup.mq20.bam

# Calculate statistics from the cleaned BAM file
samtools flagstat $name.sorted.rmdup.mq20.bam > $name.stats
```

_Realign Indels using GATK v3.1-1_

```bash
# Perform Indel Realignment
java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R CB1.fasta \
   -I $name.sorted.rmdup.mq20.bam \
   -nt 16 \
   -o $name.intervals \
   --maxIntervalSize 500 \
   --minReadsAtLocus 4

java -Xmx20g -jar GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R CB1.fasta \
   -I $name.sorted.rmdup.mq20.bam \
   -o $name.sorted.rmdup.mq20.real.bam \
   -dt NONE \
   targetIntervals $name.intervals \
   -entropy 0.075 \
   -LOD 5 \
   -maxConsensuses 30 \
   -greedy 250 \
   -model USE_READS
```
<br>

### Sorted, cleaned, realigned BAM files are ready to go!!!
