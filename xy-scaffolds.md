# Identifying putative X and Y (sex chromosome) scaffolds

_Download Bos taurus X and Y chromosome sequences_
```bash
# Download X/Y Cattle chromosomes
wget ftp://ftp.ensembl.org/pub/release-78/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna_sm.chromosome.X.fa.gz
wget ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/CHR_Y/bt_alt_Btau_4.6.1_chrY.fa.gz

# Make into Fasta file
gunzip -c Bos_taurus.UMD3.1.dna_sm.chromosome.X.fa.gz > X.fasta
gunzip -c bt_alt_Btau_4.6.1_chrY.fa.gz > Y.fasta
```

_Get the length of the X and Y, both with and without masking of repetitive elements_
```bash
# Get non-masked length of cow X
seqtk seq -A -l 0 X.fasta | grep -v "^>" | grep -o "[ATGCN]" | wc -l
	# Result: 61989626
# Get masked length of cow X
seqtk seq -A -l 0 X.fasta | grep -v "^>" | grep -o "[atgcn]" | wc -l
	# Result: 86834273

# The Y chromosome for cattle is in 12 scaffolds, concatenate them (put "N" string in between scaffolds)
cat <(echo ">Y") \
   <(seqtk seq -A -l 0 bt_alt_Btau_4.6.1_chrY.fa | \
      grep -v "^>" | \
      perl -ne 'chomp($_); print "$_"; print "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"') \
      > Y.fasta

# Get non-masked length of cow Y
seqtk seq -A -l 0 Y.fasta | grep -v "^>" | grep -o "[ATGCN]" | wc -l
	# Result: 38850781
# Get masked length of cow Y
seqtk seq -A -l 0 Y.fasta | grep -v "^>" | grep -o "[atgcn]" | wc -l
	# Result: 0
	# Result: new Y is 38,850,781 bases long
```

## X Chromosome
_Align all camel scaffolds to the cow X using LASTZ v1.02.00_
```bash
# Align using LASTZ
   # This required a cluster with at least 48GB per CPU
lastz \
	X.fasta \
	CB1.fasta \
	--step=1 \
	--gapped \
	--chain \
	--format=general \
	--verbosity=10 \
	--runtime \
	--inner=2000 \
	--ydrop=3400 \
	--gappedthresh=6000 \
	--hspthresh=2200 \
	--seed=12of19 \
	--notransition \
	--progress=10 > X.Liu.aln.tab
```

_Count number of bases in aligned scaffolds at different % coverage cutoffs_
```bash
# Loop for different cutoffs (e.g., 5 = 5% of scaffold aligned to X chromosome)
for i in 1 2 5 10 20 30 40 50 60 70 80 90
   do
   head -1 X.Liu.aln.tab >> X.aln.$i.tab
   sed '1d' X.Liu.aln.tab | \
	perl -ne 'chomp; $i='$i';@a=split(/\t/,$_); $b=$a[14]; $b=~s/%//; if ($b >= $i){print "$_\n";}' >> X.aln.$i.tab
   echo -n "Cutoff ${i}%:" >> X.cutoffs.out
   cut -f 7,9 X.aln.$i.tab | \
      sed '1d' | \
      sort | \
      uniq | \
      cut -f2 | \
      awk '{ sum+=$1} END {print sum "\t" NR}' >> X.cutoffs.out
echo "Completed $i"
done
```

_Convert all the local alignments to BED format and merge them together_
```bash
# Convert Liu LASTZ output to bed format and merge
sed '1d' X.Liu.aln.tab | \
	cut -f 7,10,11 | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i - > X.Liu.merged.bed
```

_Get the number of hits, total length, and aligned length for each camel scaffold_
```bash
# Count the number of scaffolds with alignments to the X
cut -f1 X.Liu.merged.bed | uniq > scaffolds.list
	# Results: 5158 unique scaffolds

# Build 4 column table
   # Col1 = scaffold; Col2 = total length, Col3 = # local alignments per scaffold, Col4 = aligned length
i=1
while read line
   do
   echo "Analyzing scaffold $i"
   b=$(grep "$line" CB1.fasta.fai | cut -f2)
   c=$(grep -c "$line" X.Liu.aln.tab)
   d=$(grep "$line" X.Liu.merged.bed | \
      perl -ne 'chomp;@a=split(/\t/,$_);$b=$a[2] - $a[1];print "$b\n"' | \
      awk '{ sum+=$1} END {print sum}')
   echo "$line $b $c $d" >> X.aln.scaffolds
   i=$(( $i +1 ))
done < scaffolds.list
```

_Get the mean and standard deviation of coverage per individual per scaffold_
```bash
# Make a list of the bam files (realigned and recalibrated)
ls *.sorted.rmdup.mq20.real.recal.bam > bamfiles.txt

# This code will append columns of the mean coverage for each individual bam file in the same order
# as given in the bamfiles.txt file.  Then, the standard deviation is printed in the remaining columns.
c=1
while read line
   do
   echo "Running scaffold $c"
   samtools depth -f bamfiles.txt -r $line | \
	cut -f3- | \
	Rscript mean_sd.R - | \
	tr "\n" "\t" | \
	cat - <(echo "") >> tmp.txt
   c=$(( $c + 1))
done < X.scaffolds.list

# Merge the alignments summary with coverage per individual file
paste X.aln.scaffolds tmp.txt > X.scaffolds.aln.cov
sed -i 's/ /	/g' X.scaffolds.aln.cov
rm -rf tmp.txt

# Remove scaffolds that had no depth information
cp X.scaffolds.aln.cov X.cov
cat X.cov | \
   perl -ne '$a = () = $_ =~ /\t/g; if ($a < 5){print}' | \
   cut -f1 > X.remove.list
while read line
   do
   sed -i "/$line/d" X.cov
done < X.remove.list
```

_Here is the R script 'mean_sd.R'_
```R
a = read.table("stdin")
b = colMeans(a)
c = apply(a, 2, sd)
l = dim(a)
write(as.vector(b), file = "", ncolumns = l[2], sep = "\t")
write(as.vector(c), file = "", ncolumns = l[2], sep = "\t")
```

_Using R: Calculate Wilcoxon Rank Sum test for males vs females coverage ratio_
```R
# Output a list of p-values for each scaffold
# Requires a table (file "mean.coverage") of 2 columns: col1=camel IDs col2=mean coverage
options(scipen = 999)

# Load files
a = read.table("X.cov", header = F, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
m = read.table("mean.coverage", header = F, sep = "\t")
a = t(a)

# Get subtable of scaffold coverage per individual (25 individuals)
b = a[5:29,]
b = apply(b,2,as.numeric)

# Calculate the ratio of scaffold coverage to genome-wide mean coverage
b = b/m$V2

# Set rows of known males and females
males = c(2,16,17,19,20,22,25)
females = c(9,10,11,12,13,14,15,1,18,21,23,24)

# Build my Wilcox function
my.wilcox.test <- function(x){
	test = wilcox.test(x[males], x[females], alternative = "less")
	p.val = test$p.value
	return(p.val)
}

# Apply function to each scaffold 
p.values = apply(b, 2, my.wilcox.test)
b = rbind(b, p.values)
c = rbind(as.numeric(a[2,]), as.numeric(a[3,]), as.numeric(a[4,]), b)
colnames(c) <- a[1,]
c = rbind(c, c[3,]/c[1,])
c = t(c)

# Select only scaffolds with p<0.05 & various % length cutoffs
cutoff = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
bases = vector()
for (i in cutoff){
	d = c[which(c[,29] < 0.05 & c[,30] >= i),]
	bases = c(bases, sum(d[,1]))
}

# Plot results
pdf("X.aligned.pdf", width = 7, height = 5)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(100000, 1000000000), log = "y")
box()
axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
axis(2, las = 1, at = c(100000, 1000000, 10000000, 100000000, 1000000000), labels = c(0.1, 1, 10, 100, 1000))
points(cutoff, bases, pch = 19)
points(cutoff, bases, type = "l", lwd = 1.5)
abline(h = 148823899, col = "red", lty = "dashed")
#abline(h = 61989626, col = "red", lty = "dashed") # non-repeat-masked portion
title(ylab = "Megabases", xlab = "Proportion of scaffold aligned", cex.lab = 1.4)
dev.off()

# Extract list of scaffolds for X at cutoff >20%
d = c[which(c[,29] < 0.05 & c[,30] >= 0.2),]
bases = sum(d[,1])
	# Results: Excluding 92,163,776 bases, 1,148 scaffolds
# Write list of scaffolds to exclude
write.table(names(d[,1]), file = "X.exclude", quote = F, row.names = F, col.names = F)
```

## Y Chromosome

_Align all camel scaffolds to the cow Y using LASTZ v1.02.00_
```bash
# Align using LASTZ
   # This required a cluster with at least 48GB per CPU
lastz \
	Y.fasta \
	CB1.fasta \
	--step=1 \
	--gapped \
	--chain \
	--format=general \
	--verbosity=10 \
	--runtime \
	--inner=2000 \
	--ydrop=3400 \
	--gappedthresh=6000 \
	--hspthresh=2200 \
	--seed=12of19 \
	--notransition \
	--progress=10 > Y.Liu.aln.tab
```

_Convert all the local alignments to BED format and merge them together_
```bash
# Convert LASTZ output to bed format and merge
sed '1d' Y.Liu.aln.tab | \
	cut -f 7,10,11 | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i - > Y.Liu.merged.bed
	# Results: went from 2,276,752 to 685,982 hits
```

_Get the number of hits, total length, and aligned length for each camel scaffold_
```bash
# Count the number of scaffolds with alignments to the Y
cut -f1 Y.Liu.merged.bed | uniq > Y.scaffolds.list
	# Results: 6031 unique scaffolds

# Build 4 column table
   # Col1 = scaffold; Col2 = total length, Col3 = # local alignments per scaffold, Col4 = aligned length
i=1
while read line
   do
   echo "Analyzing scaffold $i"
   b=$(grep "$line" CB1.fasta.fai | cut -f2)
   c=$(grep -c "$line" Y.Liu.aln.tab)
   d=$(grep "$line" Y.Liu.merged.bed | \
      perl -ne 'chomp;@a=split(/\t/,$_);$b=$a[2] - $a[1];print "$b\n"' | \
      awk '{ sum+=$1} END {print sum}')
   echo "$line $b $c $d" >> Y.aln.scaffolds
   i=$(( $i +1 ))
done < Y.scaffolds.list
```

_Get the mean and standard deviation of coverage per individual per scaffold_
```bash
# This will append columns of the mean coverage for each individual bam file in the same order
# as given in the bamfiles.txt file.  Then, the standard deviation is printed in the remaining columns.

c=1
while read line
   do
   echo "Running scaffold $c"
   samtools depth -f bamfiles.txt -r $line | \
	cut -f3- | \
	Rscript mean_sd.R - | \
	tr "\n" "\t" | \
	cat - <(echo "") >> tmp.txt
   c=$(( $c + 1))
done < Y.scaffolds.list

# Merge the alignments summary with coverage per individual file
paste Y.aln.scaffolds tmp.txt > Y.scaffolds.aln.cov
sed -i 's/ /       /g' Y.scaffolds.aln.cov
rm -rf tmp.txt

# Remove scaffolds that had no depth information
cp Y.scaffolds.aln.cov Y.cov
cat Y.cov | \
   perl -ne '$a = () = $_ =~ /\t/g; if ($a < 5){print}' | \
   cut -f1 > Y.remove.list
while read line
   do
   sed -i "/$line/d" Y.cov
done < Y.remove.list
	# Removed 29 scaffolds
```

_Using R: Calculate Wilcoxon Rank Sum test for males vs females coverage ratio_
```R
# Output a list of p-values for each scaffold
# Requires a table (file "mean.coverage") of 2 columns: col1=camel IDs col2=mean coverage
options(scipen = 999)

# Load files
a = read.table("Y.cov", header = F, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
m = read.table("mean.coverage", header = F, sep = "\t")
a = t(a)

# Get subtable of scaffold coverage per individual (25 individuals)
b = a[5:29,]
b = apply(b,2,as.numeric)

# Calculate the ratio of scaffold coverage to genome-wide mean coverage
b = b/m$V2

# Set rows of known males and females
males = c(2, 16, 17, 19, 20, 22, 25)
females = c(9, 10, 11, 12, 13, 14, 15, 1, 18, 21, 23, 24)

# Build my Wilcox functions; test M > F coverage
wilcox.test.mf <- function(x){
	test = wilcox.test(x[males], x[females], alternative="greater")
	p.val = test$p.value
	return(p.val)
}

# Build my Wilcox functions; test Male ratio different from 0.5
wilcox.test.m <- function(x){
	test = wilcox.test(x[males], alternative = "two.sided", mu = 0.5)
	p.val = test$p.value
	return(p.val)
}

# Build my Wilcox functions; test Female ratio less than 0.5
wilcox.test.f <- function(x){
	test = wilcox.test(x[females], alternative = "less", mu = 0.5)
	p.val = test$p.value
	return(p.val)
}

# Apply function to each scaffold 
p.values.mf = apply(b, 2, wilcox.test.mf)
	# Results: 1128 scaffolds significant
p.values.m = apply(b, 2, wilcox.test.m)
	# Results: 1144 scaffolds not different from 0.5
p.values.f = apply(b, 2, wilcox.test.f)
	# Results: 855 scaffolds significant more than 0.1
b = rbind(b, p.values.mf, p.values.m, p.values.f)
c = rbind(as.numeric(a[2,]), as.numeric(a[3,]), as.numeric(a[4,]), b)
colnames(c) <- a[1,]
c = rbind(c, c[3,]/c[1,])
c = t(c)

# Select only scaffolds with MF p<0.05 & various % length cutoffs
cutoff=c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
bases.mf = vector()
for (i in cutoff){
	d = c[which(c[,29] < 0.05 & c[,32] >= i),]
	bases.mf = c(bases.mf, sum(d[,1]))
}

# Select only scaffolds with M p>0.05 & various % length cutoffs
bases.m = vector()
for (i in cutoff){
	d = c[which(c[,30] > 0.05 & c[,32] >= i),]
	bases.m = c(bases.m, sum(d[,1]))
}

# Select only scaffolds with F p>0.05 & various % length cutoffs
bases.f = vector()
for (i in cutoff){
	d = c[which(c[,31] < 0.05 & c[,32] >= i),]
	bases.f = c(bases.f, sum(d[,1]))
}

# Select only scaffolds with F p<0.05 & M p>0.05 & various % length cutoffs
bases.final = vector()
for (i in cutoff){
	d = c[which(c[,31] < 0.025 & c[,30] > 0.025 & c[,32] >= i),]
	bases.final = c(bases.final, sum(d[,1]))
}


# Plot results
pdf("Y.aligned.pdf", height = 5, width = 7)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(100000, 1000000000), log = "y")
box()
axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
axis(2, las = 1, at = c(100000, 1000000, 10000000, 100000000, 1000000000), labels = c(0.1, 1, 10, 100, 1000))
points(cutoff, bases.final, pch = 19, col = "black")
points(cutoff, bases.final, type = "l", lwd = 1.5, col = "black")
abline(h = 38850781, col = "red", lty = "dashed")
title(ylab = "Megabases", xlab = "Proportion of scaffold aligned", cex.lab = 1.4)
dev.off()

# Extract list of scaffolds for Y at cutoff >20%
d = c[which(c[,31] < 0.025 & c[,30] > 0.025 & c[,32] >= 0.2),]
bases = sum(d[,1])
	# Results: Excluding 1,482,547 bases, 233 scaffolds
# Write list of scaffolds to exclude
write.table(names(d[,1]), file = "Y.exclude", quote = F, row.names = F, col.names = F)
```

### _Make final set of X and Y scaffolds to exclude downstream when necessary_
```bash
cat X.exclude Y.exclude > XY.exclude
   # Results:  1,381 scaffolds
```
---

### _Some extra code for making the supplmentary plots of coverage_
 _For the X_
```bash
# Get the mean depth across individuals (separate males, females, X)
# Requires a list of bamfiles and their mean coverages, females first (n=12) and males second (n=7)
# Make a new "genome" file for the X
cp CB1.fasta.fai .
c=1
while read line
   do
   echo "Grabbing scaffold $c"
   grep "$line" CB1.fasta.fai >> X.fai
   c=$(( $c + 1))
done < X.exclude

count=1
while read line
   do
   bam=$(echo -e "$line" | cut -f1)
   cov=$(echo -e "$line" | cut -f2)
   bedtools makewindows -w 1000 -s 1000 -g X.fai | \
	samtools bedcov /dev/stdin $bam | \
	perl -ne 'chomp; $c='$cov'; @a=split(/\t/,$_);$b=$a[3]/($a[2]-$a[1]);$d=$b/$c;print "$d\n"' \
	> $count.total_coverage
   count=$(( $count + 1))
done < bamfiles.sexes

paste 1.total_coverage \
	2.total_coverage \
	3.total_coverage \
	4.total_coverage \
	5.total_coverage \
	6.total_coverage \
	7.total_coverage \
	8.total_coverage \
	9.total_coverage \
	10.total_coverage \
	11.total_coverage \
	12.total_coverage \
	13.total_coverage \
	14.total_coverage \
	15.total_coverage \
	16.total_coverage \
	17.total_coverage \
	18.total_coverage \
	19.total_coverage > X.total_coverage

rm -rf {1..19}.total_coverage
```

_Graph results in R for X_
```R
# Load data
X = read.table("X.total_coverage", sep = "\t", header = F)

# Get coverage ratio across windows
means.F = rowMeans(X[,1:12])
	# Results: mean ratio of scaffold coverage to genome coverage using 1000bp non-overlapping windows: 0.9831627 (SD 0.5290527)
means.M = rowMeans(X[,13:19])
	# Results: mean ratio of scaffold coverage to genome coverage using 1000bp non-overlapping windows: 0.5708513 (SD 0.314587)

# Make plot
pdf("X.cov.pdf",width = 7, height = 5)
plot.new()
plot.window(xlim = c(0, 92753), ylim = c(0, 2))
box()
axis(1, at = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000), labels = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))
axis(2, las = 1, at = c(0, 0.5, 1, 1.5, 2), labels = c(0, 0.5, 1, 1.5, 2))
points(means.F, pch = 1, cex = 0.1, col = rgb(0, 0, 0, alpha = 0.3))
points(means.M, pch = 1, cex = 0.1, col = rgb(1, 0, 0, alpha = 0.3))
title(xlab = "Megabases", ylab = "Coverage Ratio")
dev.off()
```

_Repeat for the Y_
```bash
# Get the mean depth across individuals (separate males, females, Y)
# Requires a list of bamfiles and their mean coverages, females first (n=12) and males second (n=7)
# Make a new "genome" file for the Y
c=1
while read line
   do
   echo "Grabbing scaffold $c"
   grep "$line" CB1.fasta.fai >> Y.fai
   c=$(( $c + 1))
done < Y.exclude

count=1
while read line
do
   bam=$(echo -e "$line" | cut -f1)
   cov=$(echo -e "$line" | cut -f2)
   bedtools makewindows -w 1000 -s 1000 -g Y.fai | \
	samtools bedcov /dev/stdin $bam | \
	perl -ne 'chomp; $c='$cov'; @a=split(/\t/,$_);$b=$a[3]/($a[2]-$a[1]);$d=$b/$c;print "$d\n"' \
	> $count.total_coverage
   echo "Finished $count"
   count=$(( $count + 1))
done < bamfiles.sexes

paste 1.total_coverage \
	2.total_coverage \
	3.total_coverage \
	4.total_coverage \
	5.total_coverage \
	6.total_coverage \
	7.total_coverage \
	8.total_coverage \
	9.total_coverage \
	10.total_coverage \
	11.total_coverage \
	12.total_coverage \
	13.total_coverage \
	14.total_coverage \
	15.total_coverage \
	16.total_coverage \
	17.total_coverage \
	18.total_coverage \
	19.total_coverage > Y.total_coverage
	
rm -rf {1..19}.total_coverage
```

_Graph results in R for Y_
```R
# Load data
Y = read.table("Y.total_coverage", sep = "\t", header = F)

# Get coverage ratio across windows
means.F = rowMeans(Y[,1:12])
	# Results: mean ratio of scaffold coverage to genome coverage using 1000bp non-overlapping windows: 0.04680019 (SD 0.0908989)
means.M = rowMeans(Y[,13:19])
	# Results: mean ratio of scaffold coverage to genome coverage using 1000bp non-overlapping windows: 0.4799076 (SD 0.247982)

# Make plot
pdf("Y.cov.pdf",width = 7, height = 5)
plot.new()
plot.window(xlim = c(0, 1618), ylim = c(0, 2))
box()
axis(1, at = c(0, 250, 500, 750, 1000, 1250, 1500), labels = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5))
axis(2, las = 1, at = c(0, 0.5, 1, 1.5, 2), labels = c(0, 0.5, 1, 1.5, 2))
points(means.F, pch = 1, cex = 0.25, col = rgb(0, 0, 0, alpha = 0.5))
points(means.M, pch = 1, cex = 0.25, col = rgb(1, 0, 0, alpha = 0.5))
title(xlab = "Megabases", ylab = "Coverage Ratio")
dev.off()
```
