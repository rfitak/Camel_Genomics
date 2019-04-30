<center> <h3>Supplementary Material for:</h3> </center>
<center> <h2>Genomic signatures of domestication in Old World camels</h2> </center>

<center> <I><h5>Robert R. Fitak<sup>1</sup>, Elmira Mohandesan<sup>2</sup>, Jukka Corander<sup>3-5</sup>, Adiya Yadamsuren<sup>6</sup>, Battsetseg Chuluunbat<sup>7</sup>, Omer Abdelhadi<sup>8</sup>, Abdul Raziq<sup>9</sup>, Peter Nagy<sup>10</sup>, Chris Walzer<sup>11</sup>, Bernard Faye<sup>12</sup>, Pamela A. Burger<sup>13,14</sup></h5></I></center>

1. Department of Biology, Duke University, Durham, NC 27708 USA
2. Institute for Evolutionary Anthropology, University of Vienna, Althanstrasse 14, 1090, Vienna, Austria
3. Wellcome Trust Sanger Institute, Hinxton, UK
4. Helsinki Institute for Information Technology, Department of Mathematics and Statistics, University of Helsinki, FIN-00014, Helsinki, Finland
5. Department of Biostatistics, University of Oslo, N-0317, Oslo, Norway
6. Mammalian Ecology Laboratory, Institute of Biology, Mongolian Academy of Sciences, Peace avenue-54b, Bayanzurh district, Ulaanbaatar, 210351, Mongolia
7. Laboratory of Genetics, Institute of Biology, Mongolian Academy of Sciences, Peace avenue-54b, Bayanzurh district, Ulaanbaatar, 210351, Mongolia
8. University of Khartoum, Department for Meat Sciences, Khartoum, Sudan
9. Lasbela University of Agriculture, Water and Marine Sciences, Regional Cooperation for Development (RCD) Highway, Uthal, Pakistan
10. Farm and Veterinary Department, Emirates Industry for Camel Milk and Products, PO Box 294239, Dubai, Umm Nahad, United Arab Emirates
11. International Takhi Group - Mongolia, Baigal Ordon, Ulaanbaatar, Mongolia.
12. CIRAD-ES, UMR 112, Campus International de Baillarguet, TA C/112A, 34398, Montpellier, France
13. Research Institute of Wildlife Ecology, Vetmeduni Vienna, Savoyenstraße 1, 1160, Vienna, Austria
14. Institute of Population Genetics, Vetmeduni Vienna, Veterinaerplatz 1, 1210 Vienna, Austria

<center>
### Table of Contents
# Camel Genomics
Code for the analysis of 25 Camel genomes

[Read processing, trimming, and mapping](./read_processing.md)

Sex Chromosome Identification (TBD)

SNP Calling (TBD)

Admixture/PCA (TBD)

PSMC (TBD)

Homogeneity/HKA Tests (TBD)

Blast2GO Annotation (TBD)

[Sliding window analysis of genetic diversity](./window_diversity.md)

[GO Term Enrichment](./GO_enrichment.md)

[Recreating demographic history using G-PHOCS](./g-phocs.md) \(Currently not used\)


## Summary of the Programs/Software Used in this Study
Program | Version | Citation |
| --- | --- | --- |
| POPOOLATION | 1.2.2 | Kofler, R. et al. PoPoolation: a toolbox for population genetic analysis of next generation sequencing data from pooled individuals. PLoS One 6, e15925 (2011). |
| BWA | 0.6.2 | Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760 (2009). |
| PICARD | 1.89 | Picard Toolkit. Broad Institute, GitHub Repository. http://broadinstitute.github.io/picard/ (2018) |
| SAMTOOLS | 0.1.19 | Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079 (2009). |
| GATK | 3.1-1 | Van der Auwera, G. A. et al. From FastQ data to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Curr. Protoc. Bioinformatics 43, 11.10.11-33 (2013). |
| ANGSD | 0.563 | Korneliussen, T. S., Albrechtsen, A. & Nielsen, R. ANGSD: analysis of next generation sequencing data. BMC Bioinformatics 15, 356 (2014). |
| BEDTOOLS | <ul><li>2.17.0</li><li>2.26.0</li></ul> | Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-842 (2010). |
| VCFTOOLS | 0.1.12b | Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. 2011. The variant call format and VCFtools. Bioinformatics. 27:2156–8. |
| BEDOPS | n/a | Neph, S. et al. BEDOPS: high-performance genomic feature operations. Bioinformatics 28, 1919-1920 (2012). |
| SNPRELATE | 1.10.1 | Zheng, X. et al. A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics 28, 3326-3328 (2012). |
| ADMIXTURE | 1.23 | <ul><li>Alexander, D. H., Novembre, J. & Lange, K. Fast model-based estimation of ancestry in unrelated individuals. Genome Res. 19, 1655-1664 (2009).</li><li>Alexander, D. H. & Lange, K. Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics 12, 246 (2011).</li></ul> |
| PSMC | 0.6.4 | Li, H. & Durbin, R. Inference of human population history from individual whole-genome sequences. Nature 475, 493-496 (2011). |
| BLAST2GO | 3.0.8 | Conesa, A. et al. Blast2GO: a universal tool for annotation, visualization and analysis in functional genomics research. Bioinformatics 21, 3674-3676 (2005). |
| INTERPROSCAN | 5.7.48 | Jones, P. et al. InterProScan 5: genome-scale protein function classification. Bioinformatics 30, 1236-1240 (2014). |
| TOPGO | 2.28.0 | Alexa, A. & Rahnenfuhrer, J. (R package version 2.28.0, 2016). |
