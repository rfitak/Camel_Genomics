# Analysis of demographic history using G-phocs
Here we will use the program [G-phocs](http://compgen.cshl.edu/GPhoCS/) which was first published in [Gronau et al. 2011 *Nat Genet*](https://www.nature.com/articles/ng.937) to recreate the demographic history of the three camelid species.  G-phocs is a Bayesian coalescent sampler for inferring ancestral population sizes, population divergence times, and migration rates from individual genome sequences.  The most difficult part of the analysis is generating the input file, which must be a strict set of "neutral" loci.  Much of the G-phocs analysis below is based upon that done for gorillas ([McManus et al. 2015 *Mol Biol Evol*](https://doi.org/10.1093/molbev/msu394)) and canines ([Fan et al. 2015 *Genome Res*](http://www.genome.org/cgi/doi/10.1101/gr.197517.115))

```bash
# begin analysis code
```
