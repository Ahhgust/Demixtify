# MF files!

Demixtify prints (to standard out) a tabular file that can be parsed by
```mfPretty.R```

For example, to parse the file (named mf.tsv, made by demixtify)
```
Rscript  ~/src/Demixtify/scripts/mfPretty.R mf.tsv [ you can add multiple tsv files if you like ]
```

<br>
<br>

mfPretty.R will make a tab separated file with the following:

| Column name |  What it means              |
| ----------  |   -----------               |
| Filename    | The name of the parsed file |
| Error rate  | The mean error rate, inferred using unexpected alleles |
| NSnps       | The number of SNPs used to estimate the mixture fraction (MF) |
| MaxLL       | The maximum log likelihood. |
| MaxMF       | The MF that maximimizes the likelhood. This is the maximum likelhood estimate of the MF |
| PenultLL    | The second-highest likelihood |
| PenultMF    | The MF associated with PenultLL |
| SingleSourceLL | The likelihood associated with a MF of 0 (i.e., a single source sample) |


Note that most of the time you just care about the _MaxMF_. In the common case, the (log) likelihood ratios are so incredibly large, you needn't worry about the significance of the LLR.
<br>
<br>
Also note, the penultimate likelihoods are sometimes meaningful-- for example, if the maximum likelihood is for a single source sample (i.e., MaxMF is 0), the penultimate values are associated with a MF>0, which lets you assess a more (traditional) LR between the mixture and single source hypotheses.

## The mf.tsv file
The raw TSV produced by demixtify is a little cumbersome. It's headerless and the last row is different from the rest.

| Column number |  What it means              |
| ----------  |   -----------               |
| 1 | The label; mf means we're estimating the mixture fraction, nsnps is estimating SNP properties |
| 2 | For MF: the mixture fraction assessed ; for NSNPs, this is the estimated error rate  |
| 3 | For MF: the corresponding log (base e) likelihood ; for NSNPs, this is the number of SNPs assayed |

