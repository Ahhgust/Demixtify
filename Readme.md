# Biallelic mixture deconvolution
[Genetic genealogy](https://en.wikipedia.org/wiki/Genetic_genealogy) typically considers [SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism), often by microarray or whole genome sequencing (WGS). If your sample is a mixture, you will need a tool like *demixtify* to deconvolve your sample. The tool works by: <br>
<br>
1. **Input your mixed sample (BAM file format)
2. **Provide a B/VCF file. 
     ** This may be a "site only" VCF file (the sample is treated as having two unknown contributors)
	 ** Or it may have genotypes (and you can select which individual is a known contributor)

<br>
<br>
Demixtify will then estimate the:

*  Mixture fraction (w)
   * At a site, the fraction of reads associated with individual 1 in the mixture. Note this is *not* the same as the fraction of sampled template molecules belonging to this individual in the mixture.
   * An estimate of the genotype likelihood for each site (save those that fail QC metrics) in the VCF file.
   
Demixtify then produces a V/BCF file that is compatible with imputation software such as [GLIMPSE](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism). GLIMPSE can be used to estimate the genotype posterior probabilities. Prior to using these data it is recommended to filter/remove uncertain genotypes (e.g., with [bcftools](https://github.com/samtools/bcftools).

<br>
<br>




