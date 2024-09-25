# SNP Panels for Demixtify

Demixtify performs two analyses it:
1. * Estimates the mixture proportion
     * This is best-done using sites specifically selected for the task.
     * At present: `lowfstPanel.5bpindelGap.10maf.0.05maxfst.biallelic.sites2include.fsts.bcf`
2. * Using the estimated mixture fraction
     * It estimates joint and marginal genotypes with a fixed mixture proportion
     * This can be done on any old set of sites
     * If your goal is to use the called genotypes, you can use:
       * GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.bcf
     * If your goal is to use genotype refinement (e.g., GLIMPSE2), it is useful to consider many more sites
       * See: `hgdp1kgp_autos_sites.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.bcf`



<br>
<br>

	 


