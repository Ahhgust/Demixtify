# Biallelic mixture interpretation
When working with forensic samples, one of the first questions to consider is whether or not a sample is a mixture. <br>
Demixtify (v0.2) estimates the mixture fraction; using this mixture fraction, it then deconvolves the mixture.

### Limitations
- Demixtify is only for interpretting biallelic SNPs (as in, not indels).
  * Tri-/tetra-allelic SNPs can be split (bcftools norm -m-), but the performance of splitting needs assessment.
- Demixfity is only for deconvolving two-person mixtures.
  * Note that a three-person mixture (in truth) will better fit the hypothesis of a two-person mixture than that of a single-source sample.
- When deconvolving a very large number of sites, Demixtify can be a bit slow.
  * (Blame htslib's seek function; future versions will just stream the data)

### Recommendations

Demixtify (v0.2) is in active development. Use at your own risk.

- Either use
  * Genotype imputation (GLIMPSE or Beagle 4.1 are recommended) OR
    * EG, the output of Demixtify is the input to GLIMPSE
    * In which case, you want to type a very LARGE number of SNPs
  * Genotyping by maximum likelihood.
    * In which case considering sites in an array is recommended.
    * Filtering on the genotype quality is also recommended.
- Additionally
  * When working with "balanced" mixtures (perhaps >1:3), you will get drop-out NOT at random.
    * Really, you get drop-out depending on the genotype of the other contributor.
    * This often presents as a "matchy" profile in tools like GEDmatch. (even without imputation).
  * It should be stressed that genotype "accuracy" probably isn't the best metric to consider.



## Flags and options

Demixtify has sensible defaults. Flags/options that may (reasonably) be varied are:
```
-t nthreads (reads in the BAM file in parallel; maximum advisable value is 5)
```
<br>
<br>
```
-m min_mapping_quality (ignores reads with mapping quality < M; defaults to 20)
```
<br>
<br>

```
-b min_base_quality (ignores reads with base quality < B; defaults to 20)
```
<br>
<br>
```
-i (includes basecalls that are adjacent to an indel). Defaults to excluding such sites
```
<br>
<br>

```
-L min_read_length (ignores reads whose (mapped) length < L; default is 30)
```
<br>
<br>
Additional filters (using samflags)
```
-f read_filter (excludes reads according to SAM read filter flags). Default: 0xf04
-I read_include_filter (includes reads if all filters are met). Default: 0x2
```
<br>
See description [here](https://broadinstitute.github.io/picard/explain-flags.html) 
<br>
<br>


## Exome sequencing

We also provide files for mixture detection/deconvolution in whole exome sequencing data.<br>
See: [hg38](SupplementaryMaterial/hg38), and look for files with  "exome100bppad"  in the name.
Note, these files have only been tested for mixture detection (ie, v0.1 of the software)

<br>

## Dependencies
1. * htslib (any recent version; needed by demixtify)
     * And some (mostly) standard libraries:
       * pthreads, bzip2, zlib, libdeflate
2. * R + tidyverse (only necessary for post-processing)

## Compiling from scratch
Hopefully this won't be necessary.
And in hindsight, I used c11 threads (where I should have used c++11 threads),
which despite being part of the C11 standard, C11 threads don't have a lot of support and 
can making compiling the code difficult on older systems (ubuntu 18.04, I'm looking at you!)
*Update. I know support both c++11 threads and c11 threads. So yay for that. Older systems, as well as msys2, are now supported*
<br>
<br>
If you wish to compile the code, here is what I would suggest:
```
git clone --recursive https://github.com/Ahhgust/Demixtify.git
cd Demixtify
git submodule update --init --recursive # adds in the tokens for htslib's codecs (whatever the heck that is!)

# Let's make htslib!
#(you may be able to use a local htslib instead; just put a symlink to it in the main directory)

cd htslib 
autoreconf -i  # Build the configure script and install files it uses
./configure --disable-libcurl # somewhat lazy on my part. Libcurl can be added, but I need to move away from libhts.a to the dynamic libraries. Laziness!
make && make libhts.a # and libhts.a is the file we actually want. sooo... let's make it!
cd ..

# and let's make Demixtify!
make && echo $?

```

