# MF files!

Demixtify prints (to standard out) a tabular file that can be parsed by
```mfPretty.R```

For example, to parse the file (named mf.tsv, made by demixtify)
```
Rscript  ~/src/Demixtify/SupplementaryMaterial/scripts/mfPretty.R mf.tsv [ you can add multiple tsv files if you like ]
```

<br>
<br>

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