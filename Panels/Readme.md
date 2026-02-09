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

# hg19 support
Our analyses natively consider the hg38 reference genome. hg19 versions of (some of) the files are provided as well. In brief, the hg19 files were lifted over from hg38 (Picard); The following flags were set: <br>
--WARN_ON_MISSING_CONTIG true  --WRITE_ORIGINAL_POSITION true --WRITE_ORIGINAL_ALLELES true --RECOVER_SWAPPED_REF_ALT true  --VALIDATION_STRINGENCY LENIENT --LIFTOVER_MIN_MATCH 0.95 <br>
Additionally, the file was filtered to ensure that the hg38 and hg19 chromosomes agreed.



	 


