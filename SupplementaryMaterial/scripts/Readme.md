# Scripts

Here are some scripts that may or may not be helpful.

# vcf2panel.py
`vcf2panel.py` adds the Weir and Hill (2002) estimator of FST to to a VCF file.
It is built to recognize gnomad population labels (those named in the `bcftools annotate` command below.). <br>
An example command line is:
```
bcftools annotate -c INFO/AF_ami,AF_afr,AF_sas,AF_fin,AF_eas,AF_amr,AF_asj,AF_nfe,AN_ami,AN_afr,AN_sas,AN_fin,AN_eas,AN_amr,AN_asj,AN_nfe  -a /eva/edatums/reference_materials/gnomad/gnomad_sites_3_1_2/gnomad.genomes.v3.1.2.sites.autosomes.bcf hgdp1kgp_autos_sites.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.1maf.noindels.biallelic.sites2include.bcf | python3 ./vcf2panel.py
```
(note it write a VCF file to standard out)

