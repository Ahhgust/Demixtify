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

In brief, _Demixtify_ uses two panels; the first is a "low fst" panel. These annotations are drawn from Gnomad (v3.1); gnomad was downsampled to just SNP annotations (>1%MAF), to have a a maximum pairwise FST < 0.05, to be at least 5bp away from the nearest indel, to be intermediate in frequency (MAF > 10%), and to be accessible to whole genome sequencing (calling mask of: doi.org/10.1016/j.fsigen.2022.102785)
<br>
We *highly* recommend using this panel to infer the mixture proportion/fraction; using a different panel is not recommended without a lot of testing.
<br>
This panel was made by vcf2panel.py, which adds an estimate of FST (Weir and Hill, 2002). In our studies we found that the mean pairwise FST (i.e., average of ratios) performed better than the "overall" (i.e., ratio of averages) FST (in terms of false positives for mixture detection).

# hg19 support
Our analyses natively consider the hg38 reference genome. hg19 versions of (some of) the files are provided as well. In brief, the hg19 files were lifted over from hg38 (Picard); The following flags were set: <br>
--WARN_ON_MISSING_CONTIG true  --WRITE_ORIGINAL_POSITION true --WRITE_ORIGINAL_ALLELES true --RECOVER_SWAPPED_REF_ALT true  --VALIDATION_STRINGENCY LENIENT --LIFTOVER_MIN_MATCH 0.95 <br>
Additionally, the file was filtered to ensure that the hg38 and hg19 chromosomes agreed.


<br>

The second panel is user-defined. We provide a useful subset of the GSA (Illumina's global screening array); we recommend this for calling marginal genotypes. We also provide sites that are useful for genotype imputation/refinement. This is highly recommended if the coverage (for the relevant contributor) is low (say, <10x).

<br>

# Known limitations
Currently, Demixtify only accepts panels in the *BCF* file format. .vcf and .vcf.gz are not supported (owing to limitations of htslib).

	 


