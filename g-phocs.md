# Analysis of demographic history using G-phocs
Here we will use the program [G-phocs](http://compgen.cshl.edu/GPhoCS/) which was first published in [Gronau et al. 2011 *Nat Genet*](https://www.nature.com/articles/ng.937) to recreate the demographic history of the three camelid species.  G-phocs is a Bayesian coalescent sampler for inferring ancestral population sizes, population divergence times, and migration rates from individual genome sequences.  The most difficult part of the analysis is generating the input file, which must be a strict set of "neutral" loci.  Much of the G-phocs analysis below is based upon that done for gorillas ([McManus et al. 2015 *Mol Biol Evol*](https://doi.org/10.1093/molbev/msu394)) and canines ([Fan et al. 2015 *Genome Res*](http://www.genome.org/cgi/doi/10.1101/gr.197517.115))

1. Build a genome-wide BED file (a BED file has three columnsm 1)scaffold ID, 2) start position (0-based), and 3) end position
```bash
# Download reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Camelus_ferus/ARCHIVE/ANNOTATION_RELEASE.100/CHR_Un/cfe_ref_CB1_chrUn.fa.gz
gunzip cfe_ref_CB1_chrUn.fa.gz
mv cfe_ref_CB1_chrUn.fa CB1.fasta

# Rename scaffolds to just the accession number
sed \
   -i \
   -e "s/gi.*ref|//g" \
   -e "s/| Camelus.*//g" CB1.fasta

# Build FASTA index with SAMTOOLS
samtools faidx CB1.fasta

# Build Genome-wide BED file
bedtools \
   makewindows \
   -n 1 \
   -g CB1.fasta.fai > genome.bed
```

2.  Now we build a BED file containing regions we want to exlcude, such as repeat elements (masked regions), X and Y scaffolds, and various annotated RNAs.
```bash
# Download Annotation (Genome Release 100)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Camelus_ferus/ARCHIVE/ANNOTATION_RELEASE.100/GFF/ref_CB1_scaffolds.gff3.gz
```
