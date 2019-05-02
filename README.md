<h3><p align="center">Supplementary Methods for:</p></h3>
<h2><p align="center">Genomic signatures of domestication in Old World camels</p></h2>

<I><h5>Robert R. Fitak<sup>1,2</sup>, Elmira Mohandesan<sup>2,3</sup>, Jukka Corander<sup>4-6</sup>, Adiya Yadamsuren<sup>7,8</sup>, Battsetseg Chuluunbat<sup>9</sup>, Omer Abdelhadi<sup>10</sup>, Abdul Raziq<sup>11</sup>, Peter Nagy<sup>12</sup>, Chris Walzer<sup>13,14</sup>, Bernard Faye<sup>15</sup>, Pamela A. Burger<sup>2,14</sup></h5></I>

1. Department of Biology, Duke University, Durham, NC 27708 USA.
2. Institute of Population Genetics, Vetmeduni Vienna, Veterinärplatz 1, 1210 Vienna, Austria.
3. Institute for Evolutionary Anthropology, University of Vienna, Althanstrasse 14, 1090, Vienna, Austria.
4. Wellcome Trust Sanger Institute, Hinxton, UK.
5. Helsinki Institute for Information Technology, Department of Mathematics and Statistics, University of Helsinki, FIN-00014, Helsinki, Finland.
6. Department of Biostatistics, University of Oslo, N-0317, Oslo, Norway.
7. Institute of Remote Sensing and Digital Earth, Chinese Academy of Sciences, Jia No.20 North, DaTun road, ChaoYang District, Beijing, China.
8. Mammalian Ecology Laboratory, Institute of General and Experimental Biology, Mongolian Academy of Sciences, Peace avenue-54b, Bayanzurh District, Ulaanbaatar, 210351, Mongolia.
9. Laboratory of Genetics, Institute of General and Experimental Biology, Mongolian Academy of Sciences, Peace avenue-54b, Bayarzurh District, Ulaanbaatar, 210351, Mongolia.
 10. University of Khartoum, Department for Meat Sciences, Khartoum, Sudan.
11. Lasbela University of Agriculture, Water and Marine Sciences, Regional Cooperation for Development (RCD) Highway, Uthal, Pakistan.
12. Farm and Veterinary Department, Emirates Industry for Camel Milk and Products, PO Box 294236, Dubai, Umm Nahad, United Arab Emirates.
13. Wildlife Conservation Society, Wildlife Health Program, Bronx, New York, USA.
14. Research Institute of Wildlife Ecology, Vetmeduni Vienna, Savoyenstraße 1, 1160, Vienna, Austria.
15. CIRAD-ES, UMR 112, Campus International de Baillarguet, TA C/112A, 34398, Montpellier, France.

***
___This GitHub repository contains a summary of the various code, software, and data analysis pipelines used for the aforementioned study of 25 camel genomes. The contents represented here are only to be used as an example and not intended to be comeprehensive. The authors make no representation about the suitability or accuracy of this code, software, or data for any purpose, and make no warranties, either expressed or implied, for a particular purpose or that the use of this software or data will not infringe any third party patents, copyrights, trademarks, or other rights. The code, software and data are provided "as is". All content is hereby registered under the GNU General Public License v2.0, see [LICENSE](./LICENSE). Any publication that significantly relies upon the use of the content generated herein shall appropriately cite:___

<p align="center">Fitak et al. (in prep.) Genomic signatures of domestication in Old World camels. TBD</p>

***
  
<h2><p align="center">Table of Contents</p></h2>
<div align="center">
 
[Read processing, trimming, and mapping](./read_processing.md)

Sex Chromosome Identification (TBD)

[SNP Calling - including base and variant quality score recalibrations](./SNP-calling.md)

Admixture/PCA (TBD)

PSMC (TBD)

Homogeneity/HKA Tests (TBD)

Blast2GO Annotation (TBD)

[Sliding window analysis of genetic diversity](./window_diversity.md)

[GO Term Enrichment](./GO_enrichment.md)

[Recreating demographic history using G-PHOCS](./g-phocs.md) \(___Incomplete and not used in the study___\)

</div>

***

<h2><p align="center">Summary of the Programs/Software Used in this Study</p></h2>  

| Program | Version | Citation |
| --- | --- | --- |
| POPOOLATION | 1.2.2 | Kofler, R. et al. PoPoolation: a toolbox for population genetic analysis of next generation sequencing data from pooled individuals. PLoS One 6, e15925 (2011). |
| BWA | 0.6.2 | Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760 (2009). |
| PICARD | 1.89 | Picard Toolkit. Broad Institute, GitHub Repository. http://broadinstitute.github.io/picard/ (2018) |
| SAMTOOLS | 0.1.19 | Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079 (2009). |
| GATK | 3.1-1 | Van der Auwera, G. A. et al. From FastQ data to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Curr. Protoc. Bioinformatics 43, 11.10.11-33 (2013). |
| ANGSD | 0.563 | Korneliussen, T. S., Albrechtsen, A. & Nielsen, R. ANGSD: analysis of next generation sequencing data. BMC Bioinformatics 15, 356 (2014). |
| QUALIMAP | 2.0.1 | Okonechnikov K., Conesa A., García-Alcalde F. Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics 32(2):292-4 (2015). |
| BEDTOOLS | <ul><li>2.17.0</li><li>2.26.0</li></ul> | Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-842 (2010). |
| VCFTOOLS | 0.1.12b | Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. 2011. The variant call format and VCFtools. Bioinformatics. 27:2156–8. |
| BEDOPS | n/a | Neph, S. et al. BEDOPS: high-performance genomic feature operations. Bioinformatics 28, 1919-1920 (2012). |
| SNPRELATE | 1.10.1 | Zheng, X. et al. A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics 28, 3326-3328 (2012). |
| ADMIXTURE | 1.23 | <ul><li>Alexander, D. H., Novembre, J. & Lange, K. Fast model-based estimation of ancestry in unrelated individuals. Genome Res. 19, 1655-1664 (2009).</li><li>Alexander, D. H. & Lange, K. Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics 12, 246 (2011).</li></ul> |
| PSMC | 0.6.4 | Li, H. & Durbin, R. Inference of human population history from individual whole-genome sequences. Nature 475, 493-496 (2011). |
| BLAST2GO | 3.0.8 | Conesa, A. et al. Blast2GO: a universal tool for annotation, visualization and analysis in functional genomics research. Bioinformatics 21, 3674-3676 (2005). |
| INTERPROSCAN | 5.7.48 | Jones, P. et al. InterProScan 5: genome-scale protein function classification. Bioinformatics 30, 1236-1240 (2014). |
| TOPGO | 2.28.0 | Alexa, A. & Rahnenfuhrer, J. (R package version 2.28.0, 2016). |
