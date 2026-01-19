# MF files!

Of note, demixtify can use its own output as input (stdin). For example:

```
./demix DESCRIBE,DETECT -v lowfstPanel.5bpindelGap.10maf.0.05maxfst.biallelic.sites2include.fsts.bcf  -S - -b SyntheticMixtures/asw_asw.NA19835_NA20340.20_10.bam | ./demix DEMIX -s /dev/stdin  -v GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.bcf  -b SyntheticMixtures/asw_asw.NA19835_NA20340.20_10.bam -V Deconvolution/2unk/asw_asw.NA19835_NA20340.20_10.gsa.bcf
```

Will both estimate the MP and then deconvolve the sample based on the estimate (note the sites involved need not overlap). 


## DETECT
```
demix DETECT ...
```

Demixtify's DETECTion routine computes basic summary statistics on DESCRIBE (see below).

| Row  | What it means |
| ----------  |   ----------- |
| coverage | Average number of usable reads at the sites considered |
| snpc   | The number of (usable) sites |
| empirical_sequence_error | The mean sequence error rate, estimated (see doi.org/10.1016/j.fsigen.2023.102980) |
| reported_sequencing_error | The mean sequence error rate, reported (from base qualities) |
| mp_pointestimate  | The most likley mixture proportion, estimated by maximum likelihood |
| log_10_lr     | The log_10 likelihood ratio; 2 person vs 1 person. See also: -K/-k |
| mp_lowerbound_0.010 | lower bound on mixture fraction, Chi square approximation. See also: -a | 
| mp_upperbound_0.990 | and an upper bound ... |

Note the usage:
```
demix DESCRIBE,DETECT
```
will simply concatenate the two outputs into a single file. As describe is a 3-column file, an extra column is added to DETECT. (Simply ignore column 2)


## DESCRIBE

Demixtify estimates the mixture fraction by maximum likelihood using a grid search. Computationally, that amounts to proposing some fraction (1%), and estimating the likelihood function
given that fraction. DESCRIBE provides those correspondences, and a handful of summary statistics as well.

| Row | Value | What it means |
| --  | ----  |  ----------   |
| mf  | 0.01  | The log_e likelihood, given a MF of 0.01 |
| mf  | 0.02  | see above, but for MF=0.02. Number of steps set by: -G |
| mf  | ...   | ... |
| nsnps | Mean empirical sequence error rate | Number of SNPs evaluated |

The last line is better evaluated from `demix DETECT`
<br>
Some notes:
- You can specify the size of the grid search (`-G`); follows a uniform distribution. 
- The MF is in [0, 0.50] in the case of two unknowns
  -  As the MF approaches 0.50 (truth), the estimated MP approaches ~0.45.
  -  IE, there is some down-bias in the estimate as the MP becomes balanced.
- While it is in [0, 1.0] in the case of a known contributor.



## mfPretty

mfPretty.R is a legacy script meant to process the output from the original Demixtify (v0.1). <br>
While less pretty, consider using `demix DETECT` instead.

mfPretty.R will make a tab separated file with the following:

| Column name |  What it means              |
| ----------  |   -----------               |
| Filename    | The name of the parsed file |
| Error rate  | The mean error rate, inferred using unexpected alleles |
| NSnps       | The number of SNPs used to estimate the mixture fraction (MF) |
| Mixture_MaxLL      | The maximum log likelihood; Mixture hypothesis |
| MF          |  Mixture hypothesis; the most likely MF |
| SingleSourceLL | The likelihood associated with a MF of 0 (i.e., a single source sample) |
| LLR         | Log likelihood ratio; log_e(mixture/single source) |

If the LLR > 0 the mixture hypothesis is more likely. If it's >6.635/2 (alpha=0.01, chi square, 1 degree of freedom), the mixture hypothesis is significantly more likely. The MF provides the estimated mixture fraction.


## The mf.tsv file
The raw TSV produced by demixtify is a little cumbersome. It's headerless and the last row is different from the rest.

| Column number |  What it means              |
| ----------  |   -----------               |
| 1 | The label; mf means we're estimating the mixture fraction, nsnps is estimating SNP properties |
| 2 | For MF: the mixture fraction assessed ; for NSNPs, this is the estimated error rate  |
| 3 | For MF: the corresponding log (base e) likelihood ; for NSNPs, this is the number of SNPs assayed |

## Interpretting mixture fractions/proportions
Here's a simple framework for assessing the mixture proportion
1. If the LLR does not (significantly) favor the mixture hypothesis
   - Stop; Congratulations, treat your sample as if it were single source.
2. If the LLR *does* favor the mixture hypothesis, but the proportion is very small
   - Stop; if the proportion is small, then your sample is effectively single source (and can be treated as such). Environmental DNAs are everywhere; don't be alarmed if your sample appears to be a 0.5% mixture. In fact, the 1:199 is a special case because all the results would suggest is that the 1:199 hypothesis is more likely than the single source hypothesis. The real answer could be 1:99999!
3. If the LLR favors the mixture hypothesis and the proportion is big "enough"
   - More research needs to be done on what is big enough (and likely, it's not a question of the proportion, but also a question of the proportion AND how much data you have. IE, what the coverage of the minor contributor is), but from here you can deconvolve your sample. 