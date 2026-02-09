# Biallelic mixture deconvolution

Based on an estimated mixture proportion (taken from `-s` or `-p`), `demix DEMIX` will deconvolve the DNA mixture in the VCF file format (`-V output.vcf`).

## VCF/BCF 

### Summary statistics

Only the canonical autosomes are supported (chr1-22). The VCF header also provides:
-  the mixture proportion estimate (mixturefraction=)
-  the empirical sequence error rate (empiricalerror=)
  -  rate derived from allelic observations other than the two proposed.
-  the original bam file name (bamFile=)
-  the FST parameters (maxFstT=)
-  and the hypothesis (2unknowns or 1known+1unknown)

### Genotyping statistics

The INFO and FORMAT fields:
-  AF
  - ALT allele frequency
-  JL
  -  Joint genotype ln-likelihoods, unnormalized, in this order (Minor|Major): AA|AA,AA|AB,AA|BB,AB|AA, ... BB|BB
-  AD
  -  Allele Depth; Ref,Alt; filtered
-  OTH
  -  Number of non-ref/alt bases; no annotation means 0
-  PL
  -  Genotype likelihoods; Phred scaled, marginal *2 unknowns* or conditional *one genotype known a priori*
-  GQ
  -  Genotype Quality; Phred scaled, taken as likelihood(genotype call is right)/likelihood(wrong)

The FILTER field
-  CL
  -  Site provides conditional likelihood (1 known contributor AND they were genotyped) 
  -  Marginal likelihoods are otherwise reported


### Notes on Genotypes
Genotypes are reported for either the Major and Minor contributor (two unknowns) or for the Known/Unknown contributor (1known hypothesis). Either case provides a two-sample B/VCF file. <br>
**1known** is a bit of a misnomer, as their genotypes may be partially known (e.g., array data as known). In which case, note the `CL` tag. <br>
In addition, all likelihoods are reported, including cases where the likelihood surface is flat (e.g., 0 reads). Additional filtering on the GQ tag is recommended 
if one wishes to interpret called genotypes. Better yet, use GLIMPSE (if you need called genotypes) or use the PL tag directly.

