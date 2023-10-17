# Biallelic mixture detection
When working with forensic samples, one of the first questions to consider is whether or not a sample is a mixture. <br>
If the sample is a mixture it is also important to know the degree; balanced mixtures may need to be [deconvolved](Deconvolution.md), while if you only wish to characterize the major, you may be able to threshold the minor away.

At present, _Demixtify_ can be used to estimate the mixture fraction. Future versions will also (optionally) perform deconvolution.

<br>
<br>

1. * Input your mixed sample
     * Processed \BAM file format
       * At a minimum, duplicates should be marked (or removed).
2. * Provide a B/VCF file. 
     * This is a "sites only" VCF file
     * Some recommended files are found in [SupplementaryMaterial](SupplementaryMaterial)
       * Choose [hg38](SupplementaryMaterial/hg38) ("chr" chromosome naming, likely from UCSC) OR
       * or [grch38](SupplementaryMaterial/grch38) (no "chr"; the naming favored by ensembl)
	 

<br>
<br>
Demixtify will then estimate:

* The (log) likelihood of the mixture fraction (mf)
   * Over a range of values for mf

# Quick start (*nix systems only)
## Unix system, using static binaries


<br><br>
For examples on how to use Demixtify, see this tutorial[](SupplementaryMaterial/examples/examples.md)  

  
  
To install Demixtify, do the following:  
Clone this repository, and put it in your _src_ directory.
```
cd ~
mkdir -p src bin
cd src
git clone  https://github.com/Ahhgust/Demixtify.git
```

Put _Demixtify_ in your _bin_! (you may need to add $HOME/bin to your PATH.)
<br>
For example:

```
cp Demixtify/binaries/Nix/demix ~/bin
```

And estimate the mixture fraction
```
demix -b YOURBAM.bam -v ~/src/Demixtify/SupplementaryMaterial/hg38/GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.bcf > mf.tsv
```
_this will take a few minutes_


Demixtify also has (limited) parallel support. The computation is IO bound (ie, it takes a few seconds to estimate the MF, and a few minutes to read the bam).
You can also read the BAM file in parallel, using 4 threads for example:
```
demix -t 4 -b YOURBAM.bam -v ~/src/Demixtify/SupplementaryMaterial/hg38/GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.bcf > mf.tsv
```

Note, if you are assessing many files at once it is better to NOT use multithreading-- in the end, you're limited to the speed of the disk; once that pipe is full,
adding more threads will just make things go that much slower.

See the documentation [here](MFfile.md) to see how to work with the mixture fraction estimates.

<br>

## Flags and options

Demixtify has sensible defaults. Flags/options that may (reasonably) be varied are:<br>
For performance
```
-t nthreads (reads in the BAM file in parallel; maximum advisable value is 5)
```

Simple read filtering:
```
-m min_mapping_quality (ignores reads with mapping quality < M; defaults to 20)
```

```
-b min_base_quality (ignores reads with base quality < B; defaults to 20)
```

```
-i (includes basecalls that are adjacent to an indel). Defaults to excluding such sites
```

```
-L min_read_length (ignores reads whose (mapped) length < L; default is 30)
```

More advanced read filters (using samflags syntax)
```
-f read_filter (excludes reads according to SAM read filter flags). Default: 0xf04
-I read_include_filter (includes reads if all filters are met). Default: 0x2
```
How samflags behave is well-described by the Broad institute [here](https://broadinstitute.github.io/picard/explain-flags.html) 


Inference
```
-g grid_size (defaults to 100 (+1 for the single-source hypothesis). The grid coarseness; implicitly sets the minimum detectable fraction to 1/grid_size) 
```

```
-a AF_tag (the tag in the VCF file that provides the allele frequency; defaults to 'AF')
```



*Please ignore all other flags*! (they likely relate to deconvolution; a part of demixtify that is not ready from general use)



<br>
<br>


## Exome sequencing

We also provide files for mixture detection in whole exome sequencing data.<br>
See: [hg38](SupplementaryMaterial/hg38), and look for files with  "exome100bppad"  in the name.


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

