# Biallelic mixture detection
When working with forensic samples, one of the first questions to consider is whether or not a sample is a mixture. <br>
If the sample is a mixture it is also important to know the degree; balanced mixtures may need to be [deconvolved](Deconvolution.md), while if you only wish to characterize the major, you may be able to threshold the minor away.

At present, _Demixtify_ can be used to estimate the mixture fraction. Future versions will also (optionally) perform deconvolution.

<br>
<br>

1. * Input your mixed sample
     * Processed \BAM file format
       * At a minimum, duplicates should be marked (or removed).
2. * Provide a B/VCF file. 
     * This is a "sites only" VCF file
     * Some recommended files are found in [SupplementaryMaterial](SupplementaryMaterial)
       * Choose [hg38](SupplementaryMaterial/hg38) ("chr" chromosome naming, likely from UCSC) OR
       * or [grch38](SupplementaryMaterial/ (no "chr"; the naming favored by ensembl)
	 

<br>
<br>
Demixtify will then estimate:

* The (log) likelihood of the mixture fraction (mf)
   * Over a range of values for mf
    
<br>


<br>

## Dependencies
1. * htslib (any recent version; needed by demixtify)
2. * R + tidyverse (only necessary for post-processing)

## Exome sequencing

We also provide files for mixture detection in whole exome sequencing data.<br>
See: [hg38](SupplementaryMaterial/hg38), and look for files with  "exome100bppad"  in the name.




