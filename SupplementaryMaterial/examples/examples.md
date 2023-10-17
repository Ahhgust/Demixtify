# Examples!

It's helpful to consider some examples. We'll consider two; one that's a single source sample, and another that's a mixture.
Since these data are public, let's consider something synthetic...

# Single source data

Let's simulate some single source data!

```
bcftools view -H ../hg38/GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.bcf | python3 ../../makeRandos.py -t 4 -F AF_afr -q 30 -s 1 -r 0 | cut -f-9 | ../../demix -d /dev/stdin > singleSource4xAfr.mf
```

Which a 4x (*-t 4*) single source (*-r 0*) African/African American sample (*AF_afr*) for sites in the GSA (*bcftools view*), removing all known genotype calls (*cut*). <br>
These data are given to Demixity (*|*), which is run using the default parameters.

The format of the file we made (singleSource4xAfr.mf) is described [here](../../MFfile.md)

Looking at the first few lines of the file tells you a lot:
```
head singleSource4xAfr.mf
```

Gives a header-less TSV file (I added headers to make it more readable)


| label | mf  | log-likelihood  |
| ----- | --- | ------- |
| mf |     0.0000  |  -485744.351157961 |
| mf  |   0.0050 |  -485794.490664077 |
| mf   |   0.0100 | -485920.359570402 |
| mf   |   0.0150 | -486098.669672224 |


The basic idea of maximum likelihood techniques is to find the parameter value (mf) that maximizes the likelihood function. Log-likelihoods operate in the same way (log is a monotonic transformation, meaning that big things stay big, and small things stay small); since we're working with negative numbers, you want the log-likelihood that's closest to 0, which so happens to be the first value in the table. The associated mixture fraction is 0, suggesting that the data given is **not** a mixture. (good thing, because we know it's not!)

You can visualize the MF file ```mfPretty.R``` (note, this requires R and the tidyverse library to be installed)

![](singlesource.png =250x350)






