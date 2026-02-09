# Biallelic mixture deconvolution

Demixtify is a software suite for working with DNA mixtures. 

### Limitations
- Demixtify is only for interpretting autosomal biallelic SNPs (as in, not indels).
  * Tri-/tetra-allelic SNPs can be split (bcftools norm -m-).
- Demixtify is only for two-person mixtures.
  * Note in terms of mixture detection, the two-person hypothesis will still provide better support than the single source hypothesis (in general)


### Quick start

Grab the static binary (at present, just the x86-64 UNIX version is available),


For two unknown contributors
<br>
<br>

```
./demix -B -p LowFstPanel.BCF -v SitesToGenotype.bcf  -b YOURBAM.bam  -o Deconvolved.bcf > MFfile.tsv
```

For interpretting the MF-file, see [this](MFfile.md)
Note that the -B option uses the base-quality scores; without it, Demixtify defaults to assuming a FIXED base quality.


For one known contributor
<br>
<br>

```
./demix -B -p LowFstPanel.BCF -v SitesToGenotype.bcf  -b YOURBAM.bam  -o Deconvolved.bcf -k knownContributor.bcf > MFfile.tsv
```

The `knownContributor.bcf` is assumed to be a single-sample BCF file. If it's multisample, you can specify -K (the index of the sample; defaults to 1, the 1st individual in the BCF)


(and then use chromosome-level parallelism)

## Flags and options

Demixtify has sensible defaults. Flags/options that may (reasonably) be varied are:
```
-t nthreads (reads in the BAM file in parallel; maximum advisable value is 5)
```
<br>

```
-m min_mapping_quality (ignores reads with mapping quality < M; defaults to 20)
```
<br>

```
-b min_base_quality (ignores reads with base quality < B; defaults to 20)
```

<br>

```
-i (includes basecalls that are adjacent to an indel). Defaults to excluding such sites
```

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

## Legacy support

The original demixtify (v.10) can be found on the "legacy" branch.
See link [here](https://github.com/Ahhgust/Demixtify/tree/legacy)

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

