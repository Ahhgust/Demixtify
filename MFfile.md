# MF files!

Demixtify prints (to standard out) a tabular file that can be parsed by
```mfPretty.R```

For example, to parse the file (named mf.tsv, made by demixtify)
```
Rscript  ~/src/Demixtify/scripts/mfPretty.R mf.tsv
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


Note that most of the time you just care about the _MaxMF_. In the common case, the (log) likelihoods are so incredibly large, you needn't worry about the significance of the LLR.
<br>
<br>
Also note, the penultimate likelihoods are sometimes meaningful-- for example, if the maximum likelihood is for a single source sample (i.e., MaxMF is 0), the penultimate values are associated with a MF>0, which lets you assess a more (forensically traditional) LR between the mixture and single source hypotheses.
