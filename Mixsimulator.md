# Biallelic mixture simulator

I wrote a simple shotgun MPS DNA mixture simulator; it solely considers:

* Two person mixtures
  * And trivially, single source samples (set the mixture fraction to 0!)
* Biallelic sites
  * Though some care should be applied to sites that have been split (bcftools norm -m -)

Once you tell it the total coverage (I want a 30x genome), and the mixture fraction (9:1 mixtures, aka, a 10% mixture),
and (optionally) what population(s) you want to sample from (-E and -F flags), it will happily
generate read counts (the number of times an A allele and a B allele were observed), as well as the two sample genotypes (from the mixture)
<br>

The model is purposely very simple; it assumes a Poisson model of read depth, it only considers sequencing error, and it assumes that alleles are in  Hardy Weinberg Equilibrium.

<br>

It can, however, make output that is compatible with both VerifyBamId (pileup) and with *demixtify* (tabular).

## Quickstart

Let's start with some examples

```
bcftools view -H gnomad.1000g.phase3.100k.b38.bcf | python3 makeRandos.py -t 4 -F AF_ami -q 30 -s 1 -r 0.5 -p |VerifyBamID   --PileupFile /dev/stdin (...)
```

This makes a 4x genome (-t 4), that is a 50:50 mixture (-r 0.5) of two Amish individuals (-F AF_ami). The sequencing data are all -q 30 (0.1% sequencing error rate), and -p specifies the pileup format (needed by VerifyBamID)

```
bcftools view -H gnomad.1000g.phase3.100k.b38.bcf | python3 makeRandos.py -t 4 -F AF_ami -E AF_afr -q 30 -s 1 -r 0.8 -p |VerifyBamID   --PileupFile /dev/stdin (...)
```

This is the same as the above, except that it's a two-person mixture where the Amish person has 80% of reads on average and the other person (African/African American) has 20% of reads on average.
Note that both commands use a fixed random number seed (-s 1); this is helpful for reproducibility, however you may want to try different seeds if you are taking repeated measures.


To use with **demixtify**, instead do:

```
bcftools view -H gnomad.1000g.phase3.100k.b38.bcf | python3 makeRandos.py -t 4 -F AF_ami -E AF_afr -q 30 -s 1 -r 0.8 | cut -f-9 | demix -d /dev/stdin
```
Which treats this as a a two-person mixture with two unknowns. A single known contributor can be passed:

```
bcftools view -H gnomad.1000g.phase3.100k.b38.bcf | python3 makeRandos.py -t 4 -F AF_ami -E AF_afr -q 30 -s 1 -r 0.8 | cut -f-10 | demix -d /dev/stdin
```

which just provides one more column of information. <br><br>
```makeRandos.py``` by default makes headerless tabular data that is compatible with the BED file format (compatible with: [bedtools](https://bedtools.readthedocs.io/en/latest/) ).

| Chrom | Start Position | Stop Position | Reference Allele | Alternative Allele | Reference Counts | Alternative Counts | Other Counts | Allele frequency (reported) |  Genotype, contributor 1 | Genotype, contributor 2 |
| :---: | :---:          |         :---: | :---:            |              :---: | :---:            | :---:              | :---:        | :---:                       | :---:                    | :---:                   |
|  chr1 |  874495  | 874496 |  A  |     G   |    48 |     0  |     0  |     0.107929   |     0  |     0 |

<br>
Note that the last two columns encode genotypes as 0,1,2 (as the number of alternative alleles)
And that "other counts" refers to the number of reads supporting the other two feasible alternative alleles (these are taken as errors by *demixtify*)

