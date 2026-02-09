# SNP Panels for Demixtify

_Demixtify_ performs two analyses:
1. * Estimates the mixture proportion
     * This is best-done using sites specifically selected for the task.
     * At present: `lowfstPanel.5bpindelGap.10maf.0.05maxfst.biallelic.sites2include.fsts.bcf`
2. * Using the estimated mixture fraction
     * It estimates joint and marginal genotypes with a fixed mixture proportion
     * This can be done on any old set of sites
     * If your goal is to use the called genotypes, you can use:
       * GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.bcf
     * If your goal is to use genotype refinement (e.g., GLIMPSE2), it is useful to consider many more sites
       * See: `hgdp1kgp_autos_sites.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.chr*.bcf`

<br>

## All autosomes in 1 file.
### 30M file:
**Warning: do not use this file to DETECT or DESCRIBE!** <br>

If you want to use GLIMPSE, I highly recommend typing many many markers. Here's a link to 30758445 autosomal markers taken from the 1kGP+HGDP project.
<br>
The file is too large to post directly. To download, click on the [link](https://www.dropbox.com/scl/fi/obm0m0jguredd7b950zw2/hgdp1kgp_autos_sites.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.gnomad3_1_2anno.vcf.gz?rlkey=7ywizcysaqmqbo40ff6j0q0va&dl=1)

### 10M file
**Warning. You can use this file for DETECT or DESCRIBE, however the MP will be biased (positively)** <br>
It can still be useful, especially if the aim is to detect *balanced* mixtures in very low throughput settings (<0.20x, for example) <br>

The file is ALSO too large to post directly. To download, click on the [link](https://www.dropbox.com/scl/fi/jo81y7niosjw68t3n2bj2/hgdp1kgp_autos_sites.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.1maf.noindels.biallelic.sites2include.gnomadannos.fst.bcf?rlkey=kvw1j2e8p7j9y7iyaog4mpep8&dl=1)

# hg19 support
Our analyses natively consider the hg38 reference genome. hg19 versions of (some of) the files are provided as well. In brief, the hg19 files were lifted over from hg38 (Picard); The following flags were set: <br>
--WARN_ON_MISSING_CONTIG true  --WRITE_ORIGINAL_POSITION true --WRITE_ORIGINAL_ALLELES true --RECOVER_SWAPPED_REF_ALT true  --VALIDATION_STRINGENCY LENIENT --LIFTOVER_MIN_MATCH 0.95 <br>
Additionally, the file was filtered to ensure that the hg38 and hg19 chromosomes agreed.



	 


