#include <stdio.h>
#include <iostream>
#include <istream>
#include <string>
#include <stdlib.h>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <limits>
#include <set>
#include <chrono>

#ifdef C11THREADS
#include <threads.h> //c11 threads. which, as it turns out, are very poorly supported
#else
#include <thread> // c++11 threads, which apparently are much more standard
#endif


#include <iomanip>
#include <assert.h>
#include <random>
#include <cstdlib>

// To future travelers:
// this presumes that there is a local version of htslib in the PWD
// feel free to modify the headers below as needed
#include "htslib/htslib/sam.h"
#include "htslib/htslib/hts.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/kseq.h"

#include "demix.h"

typedef struct {
  double percentile;
  double chisq;
} ChisqTable;

// critical values from a chi-square distribution, 1 degree of freedom.
ChisqTable CHISQ[] = {
  {0.10, 2.706},
  {0.05, 5.024},
  {0.01, 6.635},
  {0.005, 7.879},
  {0.001, 10.828}
};


#define NOPE 0
#define DESCRIBE 1
#define DETECT 2
#define DEMIX 4

// sensible defaults (I hope!)
#define DEFAULT_MIN_MAPQ 20

// qualities < 20 are ignored
#define DEFAULT_MIN_BASEQ 20

// qualities > 30 are set to 30
#define DEFAULT_MAX_BASEQ 30

#define DEFAULT_MIN_READLEN 30

// size of grid search; actually is +1 (101 by default) to include the null (0% mf)
#define DEFAULT_NGRID 100

// used to exclude reads on the basis of the following:
#define DEFAULT_READ_FILTER (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)

// additionally, reads also must have ALL the following bitfields set
// (ie, must not trip the above filters, and must also be a proper pair)
#define DEFAULT_READ_INCLUDE_FILTER (BAM_FPROPER_PAIR)


// takes three likelihoods; of AA,AB,BB (note NOT phred-scaled)
// returns the most likely
#define LIKES_TO_GENO(a,b,c) (a>=b&&a>=c) ? AA : ((c>=b&&c>=a) ? BB : AB)


char summaryOutDefault[]="summary";

using namespace std;

// converts a phred quality score to its probability (no logs!)
double QUAL_2_ERRORPROB[ MAXQUAL + 1 ];

// make it so I can print a locus...
ostream& operator<<(ostream &os, const Locus &loc) {
  return os << loc.region << "\t" << loc.allele1 << "/" << loc.allele2 << "\t" << loc.af << "\t" << loc.Fst << "\t" << (unsigned)loc.genotypecall;
}


// populates the QUAL_2_ERRORPROB lookup table.
void
initQual2Prob() {
  int i;
  double *q = QUAL_2_ERRORPROB;
  for (i=0; i <= MAXQUAL; ++i, ++q)
    *q = pow(10, i/-10.0);
}




void
die(const char *arg0, const char *extra) {
	
  if (extra != NULL) {
    cerr << extra << endl;
  }
  
  cerr << "Usage " << endl <<
    arg0 << " MODE -b bamFile -v vcfFile.bcf (-k knownGenotypes.bcf) " << endl << endl <<
        "Possible modes are:" << endl << "DETECT,DESCRIBE,DEMIX" << endl << endl <<
    "Used to *detect* if a sample is a mixture" << endl <<
    "*Describe* the mixture proportion, and" << endl <<
    "*Demix* (deconvolve) the sample (VCF)" << endl << endl <<
    "Options " << endl <<

    "\t-h (prints this message) " << endl << endl <<
    "Reading/writing files:" << endl <<
    
    "\t-b bamFile or cramFile (input; Required, but see -d)" << endl <<
    "\t-T referenceFasta (input; needed to read CRAM)" << endl <<
    "\t-v/-V vcfFile (input/output; snps/snvs evaluated in the bam file)" << endl <<
    "\t-d simulation.bed (input; typically stdin; simulated WGS data; -b and -v can be omitted)" << endl <<
    "\t-s/-S summary.tsv (input (prefix) /output (filename) [from DESCRIBE]; summary statistics that characterize the mixture" << endl <<
    "\t-k knownContributor.bcf (input; a known contributor's genotypes)" << endl <<
    "\t\t-K index (input ; 1-based index of the contributor in -k)" << endl << endl <<
    "BAM processing:" << endl <<
    "\t-D r downsampling_rate (A value between [0,1]: while parsing the bam, the probability of dropping a read; defaults to 0.0, which keeps all reads)" << endl <<
    "\t-C c coverage_downsampling (A value between [0,inf]: after parsing the bam, reads are dropped to give a mean read depth of c; can combine with -D)" << endl <<
    "\t-e seed (sets the random number sEed)" << endl <<
    "\t-i (includes basecalls that are adjacent to an indel). Defaults to excluding such sites" << endl <<  
    "\t-q base_quality (minimum base quality; drops reads). Default: " << DEFAULT_MIN_BASEQ << endl <<
    "\t-Q base_quality (maximum base quality; values truncated; empirical error is bounded by this). Default: " << DEFAULT_MAX_BASEQ << endl <<
    "\t-m mapping_quality (minimum mapping quality). Default: " << DEFAULT_MIN_MAPQ <<  endl <<
    "\t-L length (minimum read length). Default: " << DEFAULT_MIN_READLEN << endl <<
    "\t-f read_filter (excludes reads according to SAM read filter flags). Default: 0x" << std::hex << DEFAULT_READ_FILTER << endl <<
    "\t-F read_include_filter (includes reads if all filters are met). Default: 0x" << std::hex << DEFAULT_READ_INCLUDE_FILTER  << endl <<
    "\t-B (tells this program to trust the quality scores; BQSR is highly recommended) " << endl <<
    "\t-R (recalibrates base qualities; BQSR of some kind is highly recommended; implies -B) " << endl <<
    "\t-t nthreads (number of threads to use. *Only used for bam/bcf (de)compression*. Defaults to 1) " << endl << 
    endl << "VCF/BCF/VCF.gz processing:" << endl <<
    "\t-W wtag (in the INFO field providing Wright's FST (Default: meanFst)" << endl<<
    "\t-A aftag (in the INFO field; the alt allele frequency. defaults to AF" << endl<<
    
    endl << "All things stats:" << endl <<
    "\t-w maxFst (discards snps with FST larger than maxFst Default: 1.0)" << endl<<
    "\t-x xmaxFst (truncates FSTs to be no larger than tmaxFst; Default: 1.0)" << endl <<
    "\t-g grid_size (for the grid-search on mixture fraction; the size of the grid (i.e., the precision)). Default: " << DEFAULT_NGRID  << endl <<
    "\t-a chisquare_index (1-5; 1: 0.10, 2: 0.05, 3: 0.01, 4: 0.005, 5: 0.001); draws confidence intervals (CIs) using pre-specified alpha values from a chi square table. Default: 3 (0.01, aka 99% CIs)" << endl <<
    "\t-p mixture_proportion (the mixture proportion; if unspecified, this is estimated)" << endl << endl <<    
 
    "For either of the filter options (-f/-F) see:" << endl <<
    "https://broadinstitute.github.io/picard/explain-flags.html" << endl << endl <<
    "Information for estimating the mixture fraction" << endl << endl;








  exit(EXIT_FAILURE);
}

bool
validNuc(char c) {
  if (c == 'A' || c== 'C' || c == 'G' || c == 'T')
    return true;
  return false;
}

bool
parseOptions(char **argv, int argc, Options &opt, char mode) {
	int i = 2;
	
	opt.filter=DEFAULT_READ_FILTER;//(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
	opt.include_filter=DEFAULT_READ_INCLUDE_FILTER;
	opt.outCounts=opt.bedFilename = opt.bamFilename= NULL;
	opt.minBaseQuality= DEFAULT_MIN_BASEQ;
	opt.maxBaseQuality= DEFAULT_MAX_BASEQ; // BQ larger than this are set to this...
	opt.minMapQuality=DEFAULT_MIN_MAPQ;
	opt.minReadLength=DEFAULT_MIN_READLEN;
	opt.outVCF="";
	opt.knowns=0;
	opt.knownBcf=NULL;
	opt.chisqIndex=2; // corresponds to alpha=0.01 
	opt.inBcf=NULL;
	opt.describeFilename=NULL;
	opt.downsampleFraction=0; // by default, keep all reads
	opt.desiredCoverage=-1; // by default, don't downsample.
	opt.referenceFasta=NULL; // only needed for CRAM support.
	opt.filterIndelAdjacent=true;
	opt.help=false;
	opt.parseVcf=false;
	opt.recalibrateQualities = opt.trustQualities=false;
	opt.maxFstD=1.; 
	opt.maxFstT=1.0; 
	
	opt.ngrid=DEFAULT_NGRID;
	opt.numThreads=1;
	opt.mixtureFraction=-1.0;
	opt.AFtag=DEFAULT_AF_TAG; // #AF
	opt.fstTag = DEFAULT_FST_TAG;
	opt.summaryFilenameRoot=string(summaryOutDefault);
	int errors=0;
	char flag;
	for ( ; i < argc; ++i) {
	  if (argv[i][0] == '-') { 
	    flag = argv[i][1];
	    if (argv[i][1] && argv[i][2]) {
	      cerr << "Sorry! Flags must be specified one at a time (as single character only). Illegal flag: " << argv[i] << endl;
	      ++errors;
	    }
	  
	  } else {
	    ++errors;
	    cerr << "Extra arguments detected." << endl;
	    break;
	  }
		
	  if (flag == 'i') {
	    opt.filterIndelAdjacent=false;
	  } else if (flag == 'h') {
	    opt.help=true;
	  } else if (flag == 'B') {
	    opt.trustQualities=true;
	  } else if (flag == 'R') {
	    opt.recalibrateQualities = opt.trustQualities=true;
	  } else { // option has an argument
	    ++i;
	    if (i==argc) {
	      cerr << "Flag " << flag << " expects and argument!" << endl;
	      ++errors;
	      break;
	    }	
	    if (flag == 'q') {
	      opt.minBaseQuality = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'Q') {
	      opt.maxBaseQuality = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'm') {
	      opt.minMapQuality = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'K') {
	      opt.knowns = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'k') {
	      opt.knownBcf = argv[i];
	      if (opt.knowns==0) // assume single-sample bcf unless otherwise specified
		opt.knowns=1;
	    } else if (flag == 'v') {
	      opt.inBcf = argv[i];
	      opt.parseVcf=true;
	    } else if (flag == 'L') {
	      opt.minReadLength = max(atoi(argv[i]), 0);
	    } else if (flag == 'f') {
	      opt.filter = atoi(argv[i]);
	    } else if (flag == 'e') {
	      int seed = atoi(argv[i]);
	      srand(seed);
	    } else if (flag == 'F') {
	      opt.include_filter = atoi(argv[i]);
	    } else if (flag == 'p') {
	      opt.mixtureFraction =   max(atof(argv[i]), 0.0000001);
	    } else if (flag == 's') {
	      opt.describeFilename = argv[i]; // output
	    } else if (flag == 'D') {
	      opt.downsampleFraction= max(atof(argv[i]), 0.);
	      if (opt.downsampleFraction>1) {
		cerr << "The downsampling fraction must be between 0 and 1!" << endl;
		++errors;
	      }
	    } else if (flag == 'C') {

	      opt.desiredCoverage=atof(argv[i]);
	      if (opt.desiredCoverage <= 0.) {
		cerr << "Coverage (-C) must be positive." << endl;
		++errors;
	      }
	      
	    } else if (flag == 'x') {
	      opt.maxFstT=max(atof(argv[i]), 0.); // maximum FST permitted ; truncates
	    } else if (flag == 'w') {
	      opt.maxFstD= max(atof(argv[i]), 0.); // maximum FST permitted ; discards markers
	    } else if (flag == 't') {
	      opt.numThreads = max(atoi(argv[i]), 1);
	    } else if (flag == 'T') {
	      opt.referenceFasta=argv[i];
	    } else if (flag == 'g') {
	      opt.ngrid = max(atoi(argv[i]), 2);
	    } else if (flag == 'V') {
	      opt.outVCF = string(argv[i]);
	    } else if (flag == 'S') {
	      opt.summaryFilenameRoot = string(argv[i]); // input
	    } else if (flag == 'b') {
	      opt.bamFilename = argv[i];
	    } else if (flag == 'A') {
	      opt.AFtag=argv[i];
	    } else if (flag == 'W') { // and the populations amongst which FST is estimated; -p implies -w
	      opt.fstTag = std::string(argv[i]);
	    } else if (flag == 'd') {
	      opt.bedFilename = argv[i]; // simulation data.
	    } else if (flag == 'a') {
	      int chisq = atoi(argv[i]);
	      if (chisq < 1 || chisq > 5) {
		cerr << "Specify alpha by indices 1-5; 1: 0.10, 2: 0.05, 3: 0.01, 4: 0.005, 5: 0.001" << endl;
		++errors;
	      } else {
		opt.chisqIndex = static_cast<unsigned char>(chisq-1);
	      }
	    } else if (flag == '-') {
		break; // respect unix-style -- to end the arguments
	    } else {
	      ++errors;
	      cerr << "Unexpected flag: " << argv[i-1] << endl;
	    }
	  }
		
	}
	
	if (opt.bamFilename==NULL) {
	  if (opt.bedFilename == NULL) { // if run using simulated data, the bed file can provide everything...
	    ++errors;
	    cerr << "No bam file was specified... I kinda need one of those!" << endl;
	  }
	}
	

	if (opt.outVCF != "" && ! (mode & DEMIX) ) {
	  cerr << "You are asking to write a VCF file (-V), but you are not selecting mode DEMIX)?!" << endl;
	  ++errors;
	}
	

	if (!opt.parseVcf && opt.bamFilename == NULL) {
	  cerr << "No BAM file provided..." << endl;
	  ++errors;
	}

	
	
	
	if (errors)
	  return false;


	if (opt.help)
	  die(argv[0], NULL);


	return true;
}

// takes a csv (str)
// and a vector (vec) of strings
// and adds the contents of str to vec
// returns the number of items added
int
csv2vec(const char* str, std::vector<std::string> *vec) {
  string tmp;
  string s(str);
  int i=0;
  stringstream ss(s);
  while (getline(ss, tmp, ',')) {
    vec->push_back(tmp);
    ++i;
  }
  return i;
}

// as above, but for doubles
int
csv2vec(const char* str, std::vector<double> *vec) {
  string tmp;
  string s(str);
  int i=0;
  stringstream ss(s);
  double d;
  char *end;
  
  while (getline(ss, tmp, ',')) {
    d = strtod(tmp.c_str(), &end);
    if (*end) {
      cerr << "Problem parsing double from: " << str << endl << end << endl;
      exit(1);
    }
    vec->push_back( d );
    ++i;
  }
  return i;
}



ostream& operator<<(ostream &os, const BaseCounter &b) {
  return os << "Ref\t" << b.refQuals.size() << "\tAlt\t" << b.altQuals.size() << "\tOther\t" << b.otherQuals.size() << "\tBad\t" << b.badCount;
}


void
initReferenceGenome(const char *fpFilename, htsFile *fp, const char *fastaRef) {

  /*
    This initializes the reference genome w/ HTSlib.
    It expects the bam/cram filename (fpFilename)
    an htsFile object (of that file)
    and the path to the reference genome (fastaRef)

    Passing this routine a BAM file is harmless, even if the fastaRef is NULL
    passing it a CRAM file will initialize things (I hope) so that the cram file can be seamlessly read 
   */
  
  int i, errors=0;

  const char *s = fpFilename;


  // this is an "ends_with" written in C
  // good times!

  for (i=0; *fpFilename; ++i, ++fpFilename)
    ;

    
  --fpFilename;
  // filenames must end in either .cram or .bam.
  // (and no, .bam is not a valid filename)
  if (i < 5) {
    cerr << "Illegal bam/cram filename: " << s << endl;
    exit(1);
  }

  if (*fpFilename != 'm')
    ++errors;
  --fpFilename;
  
  if (*fpFilename != 'a')
    ++errors;
  --fpFilename;

  if (*fpFilename == 'b') {
    return; // TODO: double check that this behavior is desired.; if it's a bam file, we don't need a reference... soo..
  } else if (! (*fpFilename == 'r' && fpFilename[-1] == 'c') )
    ++errors; // and not a cram file.
  
  if (errors) {
    cerr << "Illegal bam/cram file extension... " << s << endl;
    exit(1);
  }

  if (fastaRef==NULL) {
    cerr << "CRAM file detected, but no reference genome is provided. Please fix this!" << endl;
    exit(1);
  }

  if (hts_set_fai_filename(fp, fastaRef) != 0) {
    cerr << "Failed to load reference: " << s << endl;
    exit(1);
  }

  
}


/*
  POSSIBLE_GENOS has the genotype strings AAAA, AAAB, ... BBBB
  if we know one one of the alleles is (AB)
  this restricts POSSIBLE_GENOS to be consistent w/ what is known
  The restriction is soft; it returns the indexes in POSSIBLE_GENOS that are consistent w/ the hypothesis
  Likewise the restriction can be bsed on the known genotype benig the MAJOR, the MINOR
  EITHER (used when the minor allele fraction is estimated)
  or BOTH (trivial; nothing is known; the integers 0..8 are given
  
  Indexes are written to *itr, and the number of indexes is returned
  *itr is assumed to be of sufficient size.
  
  returns -1 on error. (type not known)

 */

int
getGenoIterators(const std::string &known, char type, int *itr) {

  int j, i=0;
  const string *genotypes=POSSIBLE_GENOS;

  // identify the indexes associated with the MAJOR contributor associated with the known genotypes
  if (type==MAJOR) {
    // see POSSIBLE_GENOS[] (demix.h) for these indexes to make sense
    // in truth, this only works with 3 possible genotypes... generalizing it is harder
    if (known[0]=='B')  // bb
      i=6;
    else if (known[1]=='B') // ab
      i=3;

    for (j=0 ; j < 3; ++j, ++itr, ++i)
      *itr = i;

    return 3;
    
  } else if (type == MINOR) {

    if (known[0]=='B')  // bb
      i=2;
    else if (known[1]=='B') // ab
      i=1;

    for (j=0 ; j < 3; ++j, ++itr, i+=3)
      *itr = i;

    return 3;

  } else if (type == EITHER) {

    // less elegant, but effective
    for ( ; i < N_GENOS; ++i, ++genotypes) {

      // major matches
      if ((*genotypes)[0] == known[0] && (*genotypes)[1] == known[1]) {
	*itr=i;
	++itr;
      } else if ((*genotypes)[2] == known[0] && (*genotypes)[3] == known[1]) { // minor matches
	*itr=i;
	++itr;
      }

    }

    return 5;
  }  else if (type==ALL) {

    for ( ; i < N_GENOS; ++i, ++itr) {
      *itr=i;
    }
    
    return N_GENOS;
  }
  return -1;
}

// converts an internal representation of a BAM sequence to ASCII; this function is just for demonstration purposes...
void
getSeq(char* buf, const uint8_t *s, int len) {
  int i;
  for (i = 0; i < len; ++i) {
    buf[i]="=ACMGRSVTWYHKDBN"[ bam_seqi(s,i) ];
  }
  buf[i] = 0;
}

/*
  Takes two coordinates (cstyle strings)
  eg chr1:100-2000
  and returns true IFF they refer to the same chromosome
 */
bool
sameChrom(const char *loc1, const char *loc2) {

  while (*loc1 && *loc2) {
    
    if (*loc1 == *loc2) {
      if (*loc1==':')
	return true;
      
    } else {
      break;
    }
    ++loc1;
    ++loc2;
  }
  return false;
}

/*
  Helper function
  region is a Locus (chr1:100-200)
  chrom is a chromosome name (chr1)
  the "match" is the chromosomes are the same.

 */
bool
matchRegToChrom(const char *reg, const char *chrom, char charEnd=':') {
  const char *p = chrom;// just a prefix match. chr1 for example
  const char *s = reg;
  while (*p && *s) {
    if (*p != *s)
      break;
    
    ++p;
    ++s;
  }
    
  if (!*p && *s==charEnd)
    return true;

  return false;
}


// sort callback
// descending order by locus index.

bool
comparePiles(const PileupReader &a, const PileupReader &b) {
  if (a.locIndex != b.locIndex)
    return b.locIndex < a.locIndex;
  
  return a.readname < b.readname;  
}



/*
  Estimates the mean read depth; an estimator of coverage.
 */
double
estimateCoverage(const BaseCounter *counts, unsigned nLoci) {
  uint64_t ntot=0;
  unsigned i;
  if (nLoci < 1) {
    cerr << "No basecalls detected." << endl;
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  for (i=0; i < nLoci; ++i, ++counts) {
    ntot += (counts->badCount + counts->otherQuals.size() + counts->refQuals.size() + counts->altQuals.size() );
  }
  return static_cast<double>(ntot)/nLoci;
}


/*
  Naive downsampling of a character vector (r)
  ds specifies the probability that a read (character) is dropped.
 */
void
downsampleVec(std::vector<unsigned char> &r, double ds) {
  unsigned i, e;
  e = r.size();
  
  for (i=0; i < r.size(); ++i) {
    if ( static_cast<double>(rand() )/RAND_MAX < ds) {
      --e;
      r[i]=r[e]; // clobber the ith data point
    }
  }
  
  r.resize( e );
    
}


/*
  Naive downsampling of a bunch of basecounter structs (vector)
  using the downsampling fraction ds (the probability that a read/character is dropped.
  downsmapling is applied to the ref and alt qualities observed, as well as the bad/other count fields.
 */
void
posthocDownsample(BaseCounter *b, unsigned nLoci,  double ds) {
  unsigned i, j, nRemove;
  
  for (j=0; j < nLoci; ++j, ++b) {

    // downsamples the ref/alt base qualities observed
    downsampleVec(b->refQuals, ds);
    downsampleVec(b->altQuals, ds);
    downsampleVec(b->otherQuals, ds);

    nRemove=0;
    for (i=0; i < b->badCount; ++i) {
      if ( static_cast<double>(rand() )/RAND_MAX < ds) 
	++nRemove;
    }
    b->badCount -= nRemove;
  }
}

/*
  Naive base recalibration. Takes in a recalibration fraction
  (ratio of means)
  reported prob a base is correct / empirical prob a base is correct
  rather than recalibrating each base quality,
  we recalibrate the lookuptable (that converts a base quality into the probability that the basecall is wrong)
 */
void
recalibrateBaseQualities(double re) {
  double *q = QUAL_2_ERRORPROB;
  for (int i=0; i <= MAXQUAL; ++i, ++q)
    *q = 1.- ((1.0- *q) * re);
  // convert error probability into correct probability, rescale it, then convert it back to an error probability
}


/*
  Converts a "pileup" into the basis of a genotype (qscores associated with ref and alt alleles; these come from Locus).
  (fills out the BaseCounter).
  Only locus locusIndex counts are evaluated.
  pile is sorted (by locus, descending, and read-pair to put read-pairs together)
  and emptied/cleared of records associated with the specified locus.
  
 */
void
processLocus(BaseCounter *count, std::vector<PileupReader> &pile, Locus *loc, uint32_t locusIndex, Options &opt) {

  uint32_t i, resizeIndex=0;
  if (pile.size() == 0)
    return;

  PileupReader *a, *b;
  char c;
  unsigned char q;
  
  count->badCount=0;
    
  // sort vector by index, then read name. paired-end overlaps will be adjacent.
  std::sort(pile.begin(), pile.end(), comparePiles);
  
  
  for (i=pile.size() -1 ; true; --i) {
    a = &pile[i];
    if (a->locIndex > locusIndex) { // sorted descending by locus; as such, implies the remainder of the records are assocaited w/ other loci
      resizeIndex=i+1;
      break;
    }
	
    if (a->locIndex != locusIndex) {
      cerr << "Assertion failed. " << a->locIndex << " " << locusIndex << endl;
      exit(1);
    }
    
    b=NULL;
    if (i>0)
      b = &pile[i-1];
      
    c = a->basecall;
    q = a->quality;
    // overlapping read-pair. 
    if (b != NULL && b->locIndex == locusIndex &&  a->readname == b->readname) {

      if (b->quality > q) { // pick the higher quality basecall
	q = b->quality;
	c = b->basecall;
      }
      --i;
    }

    // easy-peasy downsampling implemented here.
    if (opt.downsampleFraction < 1.0 &&  static_cast<double>( rand() )/RAND_MAX < opt.downsampleFraction) {
      if (i==0)
	break;
      
      continue;
    }
    if (q > MAXQUAL)
      q = MAXQUAL;
    
    if (c == loc->allele1) {
      count->refQuals.push_back( q );
    } else if (c==loc->allele2) {
      count->altQuals.push_back( q );
    } else {
      count->otherQuals.push_back( q );
    }
    
    if (i==0)
      break;
  }

  // all records are stale.
  if (resizeIndex==0) {
    pile.clear();
  } else if (resizeIndex <  pile.size()) { // last n records are stale.
    pile.resize(resizeIndex);
  }


}

uint32_t
summarizeAllRegions(samFile *in, bam1_t *b, sam_hdr_t *header, hts_idx_t *idx, std::vector<Locus> &loci, BaseCounter *counts, Options &opt) {

  int coord, unused;
  uint32_t i, currLIndex=0;
  
  hts_itr_t *iter;    // initiate an iterator
  if (loci.size() < 1) {
    cerr << "No loci sought? That doesn't sound right..." << endl;
    return 0;
  }

  for (i=0; i < loci.size(); ++i) {
    counts[i].badCount = 0;
  }
  
  std::vector<std::pair<int32_t, int64_t>> locusCoords;
  locusCoords.reserve( loci.size() );


  /* converts genomic loci (chr1:100-100) into HTSlib equivalents 
     PITA for bam, as it happens */

  int32_t chromIndex=0;
  for ( ; chromIndex < header->n_targets ; ++chromIndex) {
    bool match=  matchRegToChrom(loci[0].region.c_str(), header->target_name[ chromIndex ]);
    if (match) {
      break;
    }
  }
  
  if (chromIndex == header->n_targets) {
    cerr << "Chromosome lookup failed" << endl;
    exit(1);
  }

  
  // locusCoords are a parallel array to loci, providing the htslib equivalent to the genomic coords.
  // locusCoords are much easier to consider for operations like <
  for (std::vector<Locus>::iterator it = loci.begin(); it != loci.end(); ++it) {
    hts_parse_reg( (it->region).c_str(),  &coord, &unused);
    
    if (! matchRegToChrom((it->region).c_str(), header->target_name[ chromIndex ] )) { // false is common case.

      for ( ;! matchRegToChrom((it->region).c_str(), header->target_name[ chromIndex ] ); ++chromIndex)
	;

      if (chromIndex == header->n_targets) {
	cerr << "Chromosome lookup failed" << endl<<it->region << endl;
	exit(2);
      }   
    }

    locusCoords.emplace_back( chromIndex, coord);
  }

  
  // performance eval: this filters reads to intersect the regions listed.
  // however, is it worth the walk? Or is it better to just linearly scan?
  const char **regarray = new const char* [ loci.size() ];
  for (i=0; i < loci.size(); ++i)
    regarray[i]=loci[i].region.c_str();


  // new function; seek to all positions listed.
  iter  = sam_itr_regarray(idx, header, const_cast<char**>(regarray), loci.size() );
  if (iter==NULL) {
    cerr << "HTSlib failed to parse the regions specified" << endl << "e.g., " <<  loci[0].region << endl;
    exit(1);
  }
  
  // cerr << "Bam parsing begin" << endl;


  int offset;
  uint32_t *cigar;
  int64_t lpos, refPos, readPos, prevOp, lenCigBlock, op, len;
  
  int32_t lchromIndex; // chromosome indexes; loci must be sorted by chromIndex, position
  
  unsigned cigInd, lOffset;

  bool getNextLocus=false;

  std::vector<PileupReader> pile;
  /*
typedef struct {
  uint32_t locIndex;
  std::string readname; // deep copy of memory necessary.
  char basecall;
  unsigned char quality;
} PileupReader;
   */
  // every read that overlaps at least one locus of interest
  while ( sam_itr_next(in, iter, b) >= 0) {
    bam1_core_t core = b->core;
    
    if ( (core.flag & opt.filter) || opt.include_filter != (core.flag & opt.include_filter) ) {
      continue;
    }
    len = core.l_qseq; //length of the read.
    if (len < opt.minReadLength)
      continue;
    
    uint32_t q2 = core.qual; //mapping quality
    if (q2 < opt.minMapQuality)
      continue;

    lOffset=0;

    // locus chrom and position
    lchromIndex = locusCoords[ currLIndex].first;
    lpos = locusCoords[ currLIndex ].second;
    refPos = core.pos; // note, like coord, this is 0-based
    
        // at this point, we can only say if a locus is before the current read.
    // if so. let's clear out every locus whose information precedes this read.
    while (lchromIndex < core.tid || 
	   (core.tid==lchromIndex && lpos < refPos)) {


      processLocus( &counts[ currLIndex  ],
		    pile,
		    &(loci[currLIndex ]),
		    currLIndex,
		    opt);
      

      // one locus was consumed. No more data are possible for it.
      ++currLIndex; // double check
      if (currLIndex  >= loci.size()) {  // last locus
	break;
      }
      
      lchromIndex = locusCoords[ currLIndex].first;
      lpos = locusCoords[ currLIndex ].second;
    }
    
    if (currLIndex  >= loci.size())  // last locus
      break;

    uint8_t *s = bam_get_seq(b); //base string
    uint8_t *qual = bam_get_qual(b); // quality string (need to +33 to convert to char)

    cigar =bam_get_cigar(b);
    
    // for every locus that intersects this read...
    // common case, this doesn't loop. ie, 1 read corresponds to 1 locus.
    // since c++ doesn't support named/labeled continue statements, goto it is.
    // and for those that say goto is "bad", loops are (implemented as) gotos. So... loops are bad too?!
  nextlocus:

    if (getNextLocus) {
      
      ++lOffset;
      getNextLocus=false;
      if (currLIndex + lOffset >= loci.size()) {// last locus
	continue; // but more reads yet overlap it, so we continue
      }
    }
    
    // locus chrom and position
    lchromIndex = locusCoords[ currLIndex + lOffset ].first;
    lpos = locusCoords[ currLIndex + lOffset ].second; 
      
    prevOp=-1;

    // the read info
    refPos = core.pos;

    coord=lpos;

    for (cigInd = readPos = 0; cigInd < core.n_cigar; ++cigInd) {
      
      lenCigBlock = cigar[cigInd]>>BAM_CIGAR_SHIFT;
      op          = cigar[cigInd]&BAM_CIGAR_MASK;

      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {

	// current cigar block overlaps with the region of interest
	if (coord >= refPos && coord < (refPos + lenCigBlock) ) {

	  offset = coord - refPos;

	  // be wary of sites adjacent to alignment oddities.
	  if (opt.filterIndelAdjacent) {
	    // first base of block, treat as adjacent to indel unless there are adjacent M operators...
	    if (offset==0) {
	      if (prevOp ==-1 || !(prevOp == BAM_CMATCH ||  prevOp == BAM_CEQUAL || prevOp == BAM_CDIFF)) {

		getNextLocus=true;
		goto nextlocus;	 
	      }
	    } else if (offset == lenCigBlock-1) {

	      // last base in the block. We're skipping indel-adjacent sites; that includes if it's the last base in the (aligned portion of) the read
	      // OR 
	      // if the next cigar operation is not a "match" (including equivalent operators)
	      if (cigInd==core.n_cigar-1) {

		getNextLocus=true; // snp position overlaps, but not in a position we care to use.
		goto nextlocus;	
	      } 
	      int nextOp =  cigar[cigInd+1]&BAM_CIGAR_MASK;
	      if (nextOp != BAM_CMATCH && nextOp != BAM_CEQUAL && nextOp != BAM_CDIFF) {

		getNextLocus=true;
		goto nextlocus;	
	      }
	      
	    }
	  }
	  
	  if (offset + readPos >= len) {
	    cerr << "Problem parsing cigar string: " << core.pos << "\t" << bam_get_qname(b) << endl;
	    getNextLocus=true;
	    goto nextlocus;	
	  }
	  unsigned char qbase = qual[offset+readPos];
	  // base quality is too low...
	  if (qbase < opt.minBaseQuality) {
	    getNextLocus=true;
	    goto nextlocus;	
	  } else if (qbase > opt.maxBaseQuality)
	    qbase = opt.maxBaseQuality;

	  // not a block boundary; let's trust this value.
	  pile.emplace_back(  PileupReader{
	      (uint32_t)(currLIndex + lOffset), // which locus
		std::string( bam_get_qname(b) ), // name of read; need fresh memory alloc
		"=ACMGRSVTWYHKDBN"[bam_seqi(s,offset+readPos)], // basecall
		qbase // quality
		}
	    );

	  getNextLocus=true;
	  goto nextlocus;	
	}
		  
        refPos += lenCigBlock;
	readPos += lenCigBlock;
      } else if (op == BAM_CSOFT_CLIP) {
        readPos += lenCigBlock;
      } else if (op == BAM_CINS) {
	readPos += lenCigBlock;
      } else if (op == BAM_CDEL) {
	refPos += lenCigBlock;
      } else if (op == BAM_CREF_SKIP) {
	refPos += lenCigBlock;
      } // for our purposes, BAM_CPAD adjusts neither ther reference nor the read, ditto with hard clipping
	  
      prevOp=op;
    } // cigar parsing...


  } // next read.

  while(currLIndex  < loci.size()) {  // last locus


    processLocus( &counts[ currLIndex  ],
		  pile,
		  &(loci[currLIndex ]),
		  currLIndex,
		  opt);
      

    // one locus was consumed. No more data are possible for it.
    ++currLIndex; // double check
      
  }
    
  delete[] regarray;
  hts_itr_destroy(iter);  
  return currLIndex; // returns the number of parsed loci.
}



bool
summarizeRegion(samFile *in, bam1_t *b, sam_hdr_t *header, hts_idx_t *idx, Locus *loc, BaseCounter &count, Options &opt) {
  
  int coord, unused;
  hts_itr_t *iter;    // initiate an iterator
  iter  = sam_itr_querys(idx, header, loc->region.c_str());
  
  count.badCount=0;

  coord=unused=-1;
  if (iter==NULL) {
    cerr << "Failed to find: '" << loc->region <<"'" << endl;

    if (hts_parse_reg(loc->region.c_str(), (int*) &coord, (int*) &unused))
      fprintf(stderr, "region \"%s\" specifies an unknown reference name. Continue anyway.\n%d %d\n", loc->region.c_str(), coord, unused);
    else
      fprintf(stderr, "region \"%s\" could not be parsed. Continue anyway.\n", loc->region.c_str());

    count.badCount=SKIP_SNP;
    return false;
  }


  // coord is the 0-based coordinate start, unused is the 1-based coordinate stop
  // however we just want one base position, and 0-based
  hts_parse_reg(loc->region.c_str(),  &coord, &unused);
  
  uint32_t tot=0;  
  int offset;
  uint32_t *cigar;
  int refPos, readPos, prevOp, lenCigBlock, op, len;
  unsigned cigInd;
  
  // associate each read NAME to a base-call and quality score
  // this is so that overlapping reads can be accounted for!
  unordered_map<string, BQ> readPile;
  set<string> baddies; // in order to not double-count "bad" reads (eg, from read overlaps), let's keep a separate counting system for them

  unordered_map<string, bool> skipper;
  
  while ( sam_itr_next(in, iter, b) >= 0) {

    bam1_core_t core = b->core;
    // for paired-end/mate paired. assume mates/pairs have the same name
    string readname( bam_get_qname(b) );

    if (opt.downsampleFraction > 0) { // common case; this is 0. we're not downsampling
      
      if ( skipper.count(readname)) { // I have seen this read's mate before. the decision is already made
	if (skipper[readname]) // I skipped its mate; ergo, I skipped it.
	  continue;
	
      } else if ( static_cast<double>(rand() )/RAND_MAX < opt.downsampleFraction) { // 
	skipper[readname]=true;
	continue;
      } else
	skipper[readname]=false;
      
    }

    // none of the bad bits and all of the good bits
    if ( (core.flag & opt.filter) || opt.include_filter != (core.flag & opt.include_filter) ) {
      if (! baddies.count(readname)) {
	++count.badCount;
	baddies.insert(readname);
      }
      
      continue;
    }
    
    len = core.l_qseq; //length of the read.	
    uint8_t *s = bam_get_seq(b); //base string
    uint32_t q2 = core.qual ; //mapping quality
    uint8_t *qual = bam_get_qual(b); // quality string (need to +33 to convert to char)
    
	
    if (len < opt.minReadLength || q2 < opt.minMapQuality) {
      if (! baddies.count(readname)) {
	++count.badCount;
	baddies.insert(readname);
      }
      continue;
    }
//	getSeq(BUFF, s, len); 
//	cout << (core.flag & BAM_FREVERSE) << "\t" << core.pos << "\t" << BUFF << endl;
//	exit(1); // BAM_FREVERSE
	
    cigar =bam_get_cigar(b);
    prevOp=-1;

    for (cigInd = readPos = 0, refPos = core.pos; cigInd < core.n_cigar; ++cigInd) {
      
      lenCigBlock = cigar[cigInd]>>BAM_CIGAR_SHIFT;
      op          = cigar[cigInd]&BAM_CIGAR_MASK;

      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {

	// current cigar block overlaps with the region of interest
	if (coord >= refPos && coord < (refPos + lenCigBlock) ) {
	  
	  offset = coord - refPos;

	  // be wary of sites adjacent to alignment oddities.
	  if (opt.filterIndelAdjacent) {
	    // first base of block, treat as adjacent to indel unless there are adjacent M operators...
	    if (offset==0) {
	      if (prevOp ==-1 || !(prevOp == BAM_CMATCH ||  prevOp == BAM_CEQUAL || prevOp == BAM_CDIFF)) {
		if (! baddies.count(readname)) {
		  ++count.badCount;
		  baddies.insert(readname);
		}
		continue;
	      }
	    } else if (offset == lenCigBlock-1) {

	      // last base in the block. We're skipping indel-adjacent sites; that includes if it's the last base in the (aligned portion of) the read
	      // OR 
	      // if the next cigar operation is not a "match" (including equivalent operators)
	      if (cigInd==core.n_cigar-1) {
		if (! baddies.count(readname)) {
		  ++count.badCount;
		  baddies.insert(readname);
		}
		continue;
	      } 
	      int nextOp =  cigar[cigInd+1]&BAM_CIGAR_MASK;
	      if (nextOp != BAM_CMATCH && nextOp != BAM_CEQUAL && nextOp != BAM_CDIFF) {
		if (! baddies.count(readname)) {
		  ++count.badCount;
		  baddies.insert(readname);
		}
		continue;
	      }
	      
	    }
	  }
	  
	  if (offset + readPos >= len) {
	    cerr << "Problem parsing cigar string: " << core.pos << "\t" << bam_get_qname(b) << endl;
	    break;
	  }
	  // base quality is too low...
	  if (qual[offset+readPos] < opt.minBaseQuality) {
	    if (! baddies.count(readname)) {
	      ++count.badCount;
	      baddies.insert(readname);
	    }
	    continue;
	  }
	  
	  BQ basecall;
			  
	  basecall.b = "=ACMGRSVTWYHKDBN"[bam_seqi(s,offset+readPos)];
	  basecall.q =qual[offset+readPos];
	  if (basecall.q > MAXQUAL)
	    basecall.q = MAXQUAL;
	  
	  // same read twice (the two reads overlap)
	  // pick the best estimate of the base
	  if (readPile.count( readname)>0) {
	    if (readPile[readname].q < basecall.q) {
	      readPile[readname]=basecall;
	    }
	    
	  } else {
	    readPile[readname]=basecall;
	  }
	  break;
	}
		  
        refPos += lenCigBlock;
	readPos += lenCigBlock;
      } else if (op == BAM_CSOFT_CLIP) {
        readPos += lenCigBlock;
	//      } else if (op == BAM_CHARD_CLIP) {
      } else if (op == BAM_CINS) {
	readPos += lenCigBlock;
      } else if (op == BAM_CDEL) {
	refPos += lenCigBlock;
      } else if (op == BAM_CREF_SKIP) {
	refPos += lenCigBlock;
      } // for our purposes, BAM_CPAD adjusts neither ther reference nor the read, ditto with hard clipping
	  
      prevOp=op;
    } // cigar loop
	
    ++tot;
  }
// evaluate the pileup (correctly accounting for overlapping reads)
  for (auto iter = readPile.begin(); iter != readPile.end(); ++iter) {
    BQ basecall = iter->second;

    if (basecall.b==loc->allele1) {
      count.refQuals.push_back( basecall.q );
    } else if (basecall.b==loc->allele2) {
      count.altQuals.push_back( basecall.q );
    } else {
      count.otherQuals.push_back( basecall.q );
    }

  }

  /* performance optimization; many quality scores are the same (eg, novaseq quality binning). let's sort them
     and we can recycle computation when we estimate the genotype likelihood */
  if (count.refQuals.size() > 1)
    sort(count.refQuals.begin(), count.refQuals.end());

  if (count.altQuals.size() > 1)
    sort(count.altQuals.begin(), count.altQuals.end());
      
#if DEBUG
  cout << *loc << ' ' << count << endl;
#endif
  
  hts_itr_destroy(iter);
  return true;
}

typedef struct {
  vector<Locus> *loci;
  BaseCounter *count;
  Options opt;
  int tid;
} SummarizeRegionHelper;

int
//threadSummarizeRegion(samFile *in, bam1_t *b, sam_hdr_t *header, hts_idx_t *idx, vector<Locus> &loci, BaseCounter &count, Options &opt, int tid) {
threadSummarizeRegion(void *vargs) {
  SummarizeRegionHelper *args = (SummarizeRegionHelper*) vargs;
  int i, tid = args->tid;
  
  Options opt = args->opt;
  vector<Locus> *loci = args->loci;
  int nloc = (int) loci->size();
  BaseCounter *count = args->count;

  // borrowed from: https://dearxxj.github.io/post/5/
  bam1_t *b = bam_init1();
  samFile *in = sam_open(opt.bamFilename, "r");
  if (in == NULL)
    return -1;

  // needed for cram support...
  initReferenceGenome(opt.bamFilename, in, opt.referenceFasta);
  
  sam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL)
    return -1;
  
  hts_idx_t *idx;  // initiate a BAM index
  idx = sam_index_load(in, opt.bamFilename);   // second parameter is same as BAM file path
  if (idx == NULL) {
    cerr << "Failed to load index" << endl;
    return -1;
  }

  // statically carve up the loci
  // note number of threads <= number of loci
  // b/c the number of loci >> number of threads, static division should be optimal.
  for (i=tid; i < nloc; i += opt.numThreads) {
    summarizeRegion(in, b, header, idx, &(loci->at(i)), count[i], opt);
  }

  bam_destroy1(b);
  bam_hdr_destroy(header);
  sam_close(in);
  hts_idx_destroy(idx);
  return 1;
}


void
genotypesToAlleleWeights(double mf, double *w) {
  string *s = POSSIBLE_GENOS;
  int i=0;
  for ( ; i < N_GENOS; ++i, ++w, ++s) {
    *w = BIALLELIC_WEIGHT(*s, mf);
  }
}

// converts all proposed genotypes  e.g., (AAAB == AA individual 1, AB indivual 2)
// and the propsed mixture fraction (eg, mf==70% of reads come from individual 1)
// into the expected proportion of reads associated with an "A" allele
//
// Computes the posterior probability iff the altFrequency is defined. i.e.,  >=0)
inline void
computeLikesWithM(const BaseCounter *b, double *w, double e, double out[N_GENOS], double altFreq, double fst, bool trustLikes=false) {


  int  i = 0;
  double f;
  double like;
  const char *c;
  unsigned acount, bcount;
  
  if (!trustLikes) {
    acount = b->refQuals.size();
    bcount = b->altQuals.size();
    
    for ( ; i < N_GENOS; ++i, ++out, ++w) {
      *out = LOG_LIKE(acount, bcount, *w, e);
    }

  } else {

    for ( ; i < N_GENOS; ++i, ++out, ++w) {
      *out=0.;
      
      // formula 3 of Crysup & woerner
      // start with the bases that match the reference (bi==C in formula)
      for (const unsigned char&q : b->refQuals) {
	like = LIKE_PER_NUC(*w, q, 1);
	like += LIKE_PER_NUC(1.-*w, q, 0);
	
	*out += log(like);
      }
      // and then consider the bases that match the alternative (bi==T in formula)
      for (const unsigned char&q : b->altQuals) {
	like = LIKE_PER_NUC(*w, q, 0);
	like += LIKE_PER_NUC(1.-*w, q, 1);
	
	*out += log(like);
      }
      
      // note that the "other" class is just a constant term (and cancels out)
    }
  }
  if (altFreq >= 0.) {
    if (altFreq < MIN_AF)
      altFreq = MIN_AF;
    else if (1.-altFreq < MIN_AF)
      altFreq = 1. - MIN_AF;
	
    out -= N_GENOS;
    for (i=0 ; i < N_GENOS; ++i, ++out) {
      // todo: Fst correction?
      // and yes, given HWE I can just compute the log geno probability directly...
      c = POSSIBLE_GENOS[i].c_str();

      if (fst > 0.0)
      	f = GENOPROBFST(c, altFreq, fst);
      else
	f = GENOPROB(c, altFreq);

      c += 2;
      
      if (fst > 0.0)
	f *= GENOPROBFST(c, altFreq, fst);
      else
	f *= GENOPROB(c, altFreq);


      //      cout << POSSIBLE_GENOS[i] << " " <<  altFreq << " " << fst << "\t" << f << endl;
      *out += log(f);
    }
  }
}


// this is an estimate of the *sequencing* error rate.
// it accumulates the expected number of errors
// and divides that by the number of comparisons.
double
estimateQualError(const BaseCounter *counts, unsigned nLoci) {
  long double expectedErrors=0.;
  unsigned long ntot=0;
  for (unsigned i = 0; i < nLoci; ++i, ++counts) {
    for (auto it = (*counts).refQuals.begin(); it != (*counts).refQuals.end(); ++it) 
      expectedErrors += QUAL_2_ERRORPROB[ *it ];

    for (auto it = (*counts).altQuals.begin(); it != (*counts).altQuals.end(); ++it) 
      expectedErrors += QUAL_2_ERRORPROB[ *it ];

    for (auto it = (*counts).otherQuals.begin(); it != (*counts).otherQuals.end(); ++it) 
      expectedErrors += QUAL_2_ERRORPROB[ *it ];
    
    ntot += counts->refQuals.size() + counts->altQuals.size() + counts->otherQuals.size();
    
  }
  if (ntot < 1)
    return -1.;
  
  return expectedErrors/ntot;
  
}


// inspired by Maruki and Lynch 10.1534/g3.117.039008
// this provides a naive estimate of *sequence* error assuming a biallelic genotyping model.
// described in Crysup and Woerner 10.1534/g3.117.039008 
double
estimateError(const BaseCounter *counts, unsigned nLoci, double altMinError) {

  double e;
  unsigned i, nRight, nWrong;

  nRight=nWrong=0;  
  for (i=0; i < nLoci; ++i, ++counts) {
    nRight += counts->refQuals.size() + counts->altQuals.size();
    nWrong += counts->otherQuals.size();
  }

  if (nRight + nWrong==0) {
    return -1;
  }

  
  e = static_cast<double>(nWrong)/(nWrong+nRight);

  
  // we can only observe 2/3 possible error states
  // so lets take an expectation (assuming uniform errors)
  e *= 3.0/2;

  if (e < altMinError)
    return altMinError;
  
  return e;
}


/*
  This estimates the mixture fraction. Currently only a single-threaded version is available (
  though multithreading isn't too hard w/ this problem, it just may not save much time).

  The short of it, this does a grid search over possible mixture fractions (mf) is conducted, and
  the mf that provides the highest likelihood (assuming independence in sites and a total sampling of possible pairs of genotypes)
  is returned.

  Note that even if the mf is "known", it is known with variance (eg, you think it's a 50/50 mixture, but in truth its not)
  Also, even if it is truly a 50/50 mixture, if one of the samples is more degraded, then based on the coverage at a site
  it is not (in terms of SNPs) a 50/50 mixture.

  The MF estimator works either in the presence of 1 known genotype (which may be the major or the minor; this restricts the possible pairs of genos)
  or if both genotypes is unknown (considers all pairs of biallelic genos)

 */

// TODO: outline the likelihood estimate; keep genotypes in memory;
double
estimateMF_1Thread(const BaseCounter *counts,  vector<Locus> *loci, const Options *opt, double error, std::vector<MFest> &est, unsigned *nsnps) {

  int nLoci = (int) loci->size();
  int i,j,gridSize = opt->ngrid;
  double mf; // the proposed mixture fraction


  double loglike;
  double bestLike=-1;
  double bestMF=-1;

  double gridSizeD=gridSize;
  
  double aweights[N_GENOS];
  double genolikes[N_GENOS];
  const BaseCounter *c;
  std::vector<Locus>::iterator it; 

  int AAindexes[N_GENOS];
  int ABindexes[N_GENOS];
  int BBindexes[N_GENOS];
  int *idx;
  int nIndexes=0;
  unsigned nIncluded=0;
  bool haveKnowns = opt->knowns>0;

  
    // if the mixture fraction is specified, the the variable mf is defined
  // and the size of the grid is set to 1.
  if (opt->mixtureFraction>=0.0) {
    gridSize=0; // note that we search gridSize+1, so this ensures loops happen once.
    gridSizeD=1.;
  }
  cerr << "estimate MF?\n";
  
  if (haveKnowns) {
    // nIndexes==3 in all cases.
    nIndexes = getGenoIterators("AA", MINOR, AAindexes);
    getGenoIterators("AB", MINOR, ABindexes);
    getGenoIterators("BB", MINOR, BBindexes);

  } else { // both unknowns, then the mixture fraction only meaningfully varies from 0-0.5 (or from 0.5-1, take you pick!)
    gridSizeD += gridSizeD;
  }
  
  for (i = 0; i <= gridSize; ++i) {
    loglike=0.;
    
    if (opt->mixtureFraction>= 0.0)  // happens once; breaks at end of loop
      mf = opt->mixtureFraction;
    else // actual looping...
      mf = i/gridSizeD;

    // convert the 9 genotype calls with the proposed mixture fraction
    // to get the fraction of reads expected of an "A" allele (a continuous "genotype")
    genotypesToAlleleWeights( mf, aweights);
    
    it = loci->begin();
    c = counts;
    for (j=0; j < nLoci; ++j, ++it, ++c) {

      if (c->badCount < IGNORE_SNP && c->refQuals.size() + c->altQuals.size() > 0) {

	computeLikesWithM(c, aweights, error, genolikes, it->af, it->Fst,opt->trustQualities);
	
	if (i==0)
	  ++nIncluded;
	
	if (haveKnowns) {
	  if (it->genotypecall==NOCALL)
	    continue;
	  else if (it->genotypecall==AA) 
	    idx=AAindexes;
	  else if (it->genotypecall==AB) 
	    idx=ABindexes;
	  else if (it->genotypecall==BB) 
	    idx=BBindexes;
	  else {
	    cerr << "Should never happen! " << *it << endl;
	    exit(1);
	  }

	  loglike += log_sum_exp_with_knowns(genolikes, idx, nIndexes);

	} else {
	  double ll = log_sum_exp(genolikes);
	  /*
	  if (ll > 0) {
	    cerr << ll << " "  << c->refCount << " " <<  c->altCount << "\t" << it->af << endl;
	    for (int x =0; x <9; ++x)
	      cerr << "\t" << genolikes[x];
	  }
	  cerr << endl;
	  */

	  loglike += ll;
	}
	
      }
    }

    if (i==0 || bestLike < loglike) {
      bestLike=loglike;
      bestMF=mf;
    }
    //cout << mf << "\t" << loglike << endl;
    // printf("mf\t%0.4f\t%0.9f\n", mf, loglike);
    est[i].mf=mf;
    est[i].loglike=loglike;
    if (opt->mixtureFraction>= 0.0)  // happens once; breaks at end of loop
      break;
  }

  
  //printf("nsnps\t%f\t%u\n", error, nIncluded);
  *nsnps=nIncluded;
  return bestMF;
}

// logsumexp trick
// applied to the genos;
// used to get the log likelihood of P(D|mf) = product over SNP( sum over genotypes ( a.count,c.count|mf, g))
inline double
log_sum_exp(double *genolikes) {
  int i;
  double s=0.;
  double max = *genolikes;
  double *g = genolikes+1;


  for (i=1; i < N_GENOS; ++i, ++g) {
    if (*g > max)
      max = *g;

  }
  
  for (i=0; i < N_GENOS; ++i, ++genolikes) 
    s += exp(*genolikes-max);

  
  return log(s) + max;
  
}

// same as log_sum_exp
// however only the indexes in *knowns are traversed (all ngenos of them)
inline double
log_sum_exp_with_knowns(double *genolikes, int *knowns, int ngenos) {
  int i, *k=knowns+1;
  double s=0.;
  double max = genolikes[ *knowns ];


  for (i=1; i < ngenos; ++i, ++k) {
    if (genolikes[ *k ] > max)
      max = genolikes[ *k ];
  }
  
  k=knowns;

  for (i=0; i < ngenos; ++i, ++k) 
    s += exp(genolikes[*k]-max);

  
  return log(s) + max;
  
}


/*
  converts genotype log-likelihoods (ln(likelihood))
  into Phred-scaled (and normalized) likelihoods (-10*log10(likelihood))-(-10*log10(max(likelihood)))
  see: https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
  and
  https://en.wikipedia.org/wiki/Phred_quality_score
  ie, this will make the values that go into the PL tag in the VCF file format.
 */
void
logLikeToPL(const double *in, int n, double *out) {

  int i;
  // the maximum likelhood estimate; phred-scaled
  double maxLikePl = (*std::max_element(in, in +n)) * PHRED_SCALAR;
  
  for(i=0; i < n; ++i, ++in, ++out)
    *out = (*in)*PHRED_SCALAR - maxLikePl;
  
}

void
logLikeToPL(const double *in, int n, int *out) {

  int i;
  // the maximum likelhood estimate; phred-scaled
  double maxLikePl = (*std::max_element(in, in +n)) * PHRED_SCALAR;
  
  for(i=0; i < n; ++i, ++in, ++out)
    *out = round((*in)*PHRED_SCALAR - maxLikePl);
  
}

/*
  Syntactic sugar.
  This reduces the 9 joint genotypes likelihoods all pairs of : (AA|AB|BB)
  into the 6 marginal genotypes; (AA|AB|BB) for the major and the minor
  *out must be of size 6.
 */
void
computeMarginals(int AAindexes[N_GENOS],
	     int ABindexes[N_GENOS],
	     int BBindexes[N_GENOS],
	     int AAindexesMinor[N_GENOS],
	     int ABindexesMinor[N_GENOS],
	     int BBindexesMinor[N_GENOS],
	     double *genolikes,
	     double *out) {

  
  *out++ = log_sum_exp_with_knowns(genolikes, AAindexesMinor, 3);
  *out++ = log_sum_exp_with_knowns(genolikes, ABindexesMinor, 3);
  *out++ = log_sum_exp_with_knowns(genolikes, BBindexesMinor, 3);
  *out++ = log_sum_exp_with_knowns(genolikes, AAindexes, 3);
  *out++ = log_sum_exp_with_knowns(genolikes, ABindexes, 3);
  *out = log_sum_exp_with_knowns(genolikes, BBindexes, 3);


}

//
// Takes a vector of length 3; log likelihoods for AA,AB,BB genotypes
// and returns the Phred-scaled genotype quality (int)
// min/max set to [0,99], ala GATK
int
getGQ(double *marginalLikes) {

  int maxI=0; //argmax
  int i;
  double probWrong=0.;
  for (i=1; i < 3; ++i) {
    if (marginalLikes[i] > marginalLikes[maxI])
      maxI=i;
  }

  // gets the sum of the likelihoods associated w/ all other genotype calls considered
  for(i=0; i < 3; ++i) {
    if (i != maxI)
      probWrong += exp(marginalLikes[i]);
  }

  //cout << exp(marginalLikes[maxI]) << "\t" << probWrong << "\t" << PHRED_SCALAR*(log(probWrong)-marginalLikes[maxI]) << endl;
  
  return round(std::min(99., std::max(0., (PHRED_SCALAR*(log(probWrong)-marginalLikes[maxI])))));
}


unsigned 
populateLoci(vector<Locus> *loci, Options &opt, string chromosome) {


  cerr << opt.bedFilename  << " " << chromosome << endl;
  VcfIterator wgsPanel = initVcfIterator(opt.bedFilename, 0, true);
  hts_itr_t *itr = bcf_itr_querys(wgsPanel.idx, wgsPanel.hdr, chromosome.c_str());

  int ret;
  loci->clear();

  while ((ret=bcf_itr_next(wgsPanel.fp, itr, wgsPanel.rec)) == 0) {

    bcf_unpack(wgsPanel.rec, BCF_UN_STR);

    if (!bcf_is_snp(wgsPanel.rec) || wgsPanel.rec->n_allele > 2)
      continue;

    
    string ch( bcf_hdr_id2name(wgsPanel.hdr, wgsPanel.rec->rid));

    if (ch != chromosome)
      break;
    
    Locus loc;
    
    // make a locus object
    loc.region=ch + ":" + std::to_string(wgsPanel.rec->pos+1) + "-" + std::to_string(wgsPanel.rec->pos+1);
    loc.allele1 = wgsPanel.rec->d.allele[0][0];
    loc.allele2 = wgsPanel.rec->d.allele[1][0];
    loc.pos = wgsPanel.rec->pos+1;
    loc.genotypecall=NOCALL;

    // can consider emplace_back if this is too slow.
    loci->push_back(loc);
  }
  cerr << loci->size() << " " << ret << endl;
  
  closeVcfIterator(wgsPanel);
  return loci->size();
}


void
closeVcfWriter(VcfWriter &writer, const Options &opt) {
  if (writer.fp!= NULL) {
    if (hts_close(writer.fp))
      cerr << "Failed to close " << writer.filename << endl;

  }

  if (writer.hdr !=NULL) {
    bcf_hdr_destroy(writer.hdr);

  } else
    return;

  if (writer.file_format == "wb") {
    // index
    if ( bcf_index_build3(writer.filename.c_str(), NULL,14,opt.numThreads))
      cerr << "Some problem writing the B/VCF file index. That's not good!" << endl;
  }
}

// heavily influenced by: https://github.com/odelaneau/GLIMPSE/blob/master/phase/src/io/genotype_writer.cpp
// see also: http://broadinstitute.github.io/gamgee/doxygen/hts_8h.html
void
writeVcfOneChrom(double mf, double error, vector<Locus> *loci, Options &opt, const std::string &chrom, VcfWriter *writer) {


  unsigned len, i;

  // first time through, the writer struct is initialized
  // and a generic autosomal header is written
  if (writer->fp==NULL) {
    writer->file_format = "wb";
    writer->filename=opt.outVCF;
    
    if (opt.outVCF != "-") {
      len = opt.outVCF.length();
      if (len > 4 && opt.outVCF[len-1] == 'f' &&  opt.outVCF[len-2] == 'c'  && opt.outVCF[len-3] == 'v')
	writer->file_format = "w"; // uncompressed VCF if that's what's asked for.
    }
	      
    writer->fp = hts_open(opt.outVCF.c_str(),writer->file_format.c_str());
    
    if (opt.numThreads > 1)
      hts_set_threads(writer->fp, opt.numThreads);


    writer->hdr = bcf_hdr_init("w"); 
    
    bcf_hdr_append(writer->hdr, string("##source=Demixtify v" + std::string(VERSION)).c_str());
    bcf_hdr_append(writer->hdr, string("##mixturefraction=" + std::to_string(mf)).c_str()); // yay c++11!
    bcf_hdr_append(writer->hdr, string("##empiricalerror=" + std::to_string(error)).c_str());
    // log key items from how this file was made
    bcf_hdr_append(writer->hdr, string("##bamFile=" + std::string(opt.bamFilename)).c_str() );
    bcf_hdr_append(writer->hdr, string("##maxFstT=" + std::string(opt.fstTag) + ";" +  std::to_string(opt.maxFstT)).c_str() );

  // and clearly state hypotheses (along w/ reminders as to what was used)
    if (opt.knowns) {
      bcf_hdr_append(writer->hdr, string("##hypothesis=1known+1unknown").c_str() );
      bcf_hdr_append(writer->hdr, string("##knownIndex=" + std::to_string(opt.knowns)).c_str() );
      bcf_hdr_append(writer->hdr, string("##knownFile=" + std::string(opt.knownBcf)).c_str() );
    } else
      bcf_hdr_append(writer->hdr, string("##hypothesis=2unknowns").c_str() );


    // TODO; expand to X?
    for (i=1; i < 23; ++i) {
      std::string thischrom="chr" + to_string(i);
	
      bcf_hdr_append(writer->hdr, std::string("##contig=<ID="+ thischrom  + ">").c_str());

    }

    bcf_hdr_append(writer->hdr, "##FILTER=<ID=CL,Description=\"Site provides conditional likelihood\">");
    bcf_hdr_append(writer->hdr, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"ALT allele frequency\">");
    bcf_hdr_append(writer->hdr, "##INFO=<ID=JL,Number=9,Type=Float,Description=\"Joint genotype ln-likelihoods, unnormalized, in this order (Minor|Major): AA|AA,AA|AB,AA|BB,AB|AA, ... BB|BB, \">");  
    bcf_hdr_append(writer->hdr, "##INFO=<ID=AD,Number=2,Type=Integer,Description=\"Allele Depth; Ref,Alt; filtered\">");
    bcf_hdr_append(writer->hdr, "##INFO=<ID=OTH,Number=1,Type=Integer,Description=\"Number of non-ref/alt bases; no annotation means 0\">");
  
    bcf_hdr_append(writer->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotypes\">");
    bcf_hdr_append(writer->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods; Phred scaled, marginal\">");
    bcf_hdr_append(writer->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality; Phred scaled, taken as likelihood(genotype call is right)/likelihood(wrong)\">");


    if (!opt.knowns) {
      bcf_hdr_add_sample(writer->hdr, "Major");
      bcf_hdr_add_sample(writer->hdr, "Minor");
    } else  if (mf > 0.5) { // minor is unknown
      bcf_hdr_add_sample(writer->hdr, "Unknown");
      bcf_hdr_add_sample(writer->hdr, "Known");

    } else {
      bcf_hdr_add_sample(writer->hdr, "Known");
      bcf_hdr_add_sample(writer->hdr, "Unknown");
    }
    
    bcf_hdr_add_sample(writer->hdr, NULL);      // to update internal structures ; may or may not be necessary?!
  
    if (bcf_hdr_write(writer->fp, writer->hdr) ) {
      cerr << "Failed to write header to file: " << opt.outVCF << endl;
      exit(1);
    }
  }

  bcf1_t *rec = bcf_init1();

  /* data structures for converting allele counts into mixed/unmixed genotypes */
  int AAindexes[N_GENOS];
  int ABindexes[N_GENOS];
  int BBindexes[N_GENOS];

  int AAindexesMinor[N_GENOS];
  int ABindexesMinor[N_GENOS];
  int BBindexesMinor[N_GENOS];
  
  int *idx;

  double marginalLikes[N_MARGINALS];
  
  getGenoIterators("AA", MAJOR, AAindexes);
  getGenoIterators("AB", MAJOR, ABindexes);
  getGenoIterators("BB", MAJOR, BBindexes);
  
  getGenoIterators("AA", MINOR, AAindexesMinor);
  getGenoIterators("AB", MINOR, ABindexesMinor);
  getGenoIterators("BB", MINOR, BBindexesMinor);

  // convert genotypes + mixture fractions to allele weights
  double aweights[N_GENOS];
  genotypesToAlleleWeights( mf, aweights);

  double genolikes[N_GENOS];
  float genolikesf[N_GENOS];
  int phredlikes[N_MARGINALS];


  //Add records
  // 5: 2 diploid samples (4) + 1 null record (bcf_int32_vector_end)
  int bcfgenotypes[5];
  int clID = bcf_hdr_id2int(writer->hdr,  BCF_DT_ID, "CL");

  BaseCounter counts;

  // Stuff we need to process bam 
  bam1_t *b = bam_init1();
  samFile *in = sam_open(opt.bamFilename, "r");
  if (in == NULL)
    exit(EXIT_FAILURE);

  if (opt.numThreads>1)
    hts_set_threads(in, opt.numThreads);
  
  initReferenceGenome(opt.bamFilename, in, opt.referenceFasta);
      
  sam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL)
    exit(EXIT_FAILURE);
  
  hts_idx_t *bamidx;  // initiate a BAM index
  bamidx = sam_index_load(in, opt.bamFilename);   // second parameter is same as BAM file path
  if (bamidx == NULL) {
    cerr << "Failed to load index for " << opt.bamFilename << endl;
    exit(EXIT_FAILURE);
  }

  uint32_t chromIndex = bcf_hdr_name2id(writer->hdr, chrom.c_str() );   

  for (std::vector<Locus>::iterator it = loci->begin(); it != loci->end(); ++it) {

    // clear datums...
    bcf_clear1(rec);
    counts.refQuals.clear();
    counts.altQuals.clear();
    counts.otherQuals.clear();
    counts.badCount=0;

    
    rec->rid = chromIndex;
    rec->pos = it->pos-1; // convert from 1- to 0-based coordinate
    
    // add ref and alt allele categories
    std::string alleles = std::string(1, it->allele1) + "," + std::string(1, it->allele2);
    bcf_update_alleles_str(writer->hdr, rec, alleles.c_str());
    
    
    // fill in the counts data (bam processing)
    if (! summarizeRegion(in, b, header, bamidx, &(*it), counts, opt)) {
      cerr << "Skipping region " << it->region << endl; 
      continue;
    }
    
#if DEBUG
    cerr << *it << endl;
    cerr << counts << endl;
#endif
    
    // default genotypes set to 0/0 (common case)
    std::fill(bcfgenotypes, bcfgenotypes + 4, bcf_gt_unphased(false));
    bcfgenotypes[4] = bcf_int32_vector_end;
    
    // -1s turn off posterior prob calc
    // todo: use equation 3
    computeLikesWithM(&counts, aweights, error, genolikes, -1, -1, opt.trustQualities);

    if (!opt.knowns || it->genotypecall==NOCALL ) {
      
      // corner case; minor is known, in which case the joint genotype matrix is flipped
      if (mf > 0.5)
	computeMarginals(AAindexesMinor,ABindexesMinor,BBindexesMinor,AAindexes,ABindexes,BBindexes,genolikes,marginalLikes);
      else
	computeMarginals(AAindexes,ABindexes,BBindexes,AAindexesMinor,ABindexesMinor,BBindexesMinor,genolikes,marginalLikes);
      
      //computeMarginals(AAindexes,ABindexes,BBindexes,AAindexesMinor,ABindexesMinor,BBindexesMinor,genolikes,marginalLikes);
      //logLikeToPL(const double *in, int n, int *out) {
      logLikeToPL(marginalLikes, 3, phredlikes); 
      logLikeToPL(marginalLikes+3, 3, phredlikes+3); 

      int g;
      g = LIKES_TO_GENO(marginalLikes[0], marginalLikes[1], marginalLikes[2]);
	
      if (g == AB) {
	bcfgenotypes[1] = bcf_gt_unphased(true);
      } else if (g==BB) {
	bcfgenotypes[0] = bcf_gt_unphased(true);
	bcfgenotypes[1] = bcf_gt_unphased(true);
      }

      g = LIKES_TO_GENO(marginalLikes[3], marginalLikes[4], marginalLikes[5]);      
      
      if (g == AB) {
	bcfgenotypes[3] = bcf_gt_unphased(true);
      } else if (g==BB) {
	bcfgenotypes[2] = bcf_gt_unphased(true);
	bcfgenotypes[3] = bcf_gt_unphased(true);
      }

      
    } else {

      bcf_update_filter(writer->hdr, rec, &clID, 1);
      // and given what we know, grab the appropriate set of 3 genotypes of the joint likelihoods,

      idx = AAindexesMinor;
      if (it->genotypecall == AB) {
	idx = ABindexesMinor;
      } else if (it->genotypecall == BB) {
	idx = BBindexesMinor;
      }
      

      for (i=0; i < 3; ++i, ++idx) {
	marginalLikes[i] = genolikes[ *idx ];
      }
      
      int g = LIKES_TO_GENO(marginalLikes[0], marginalLikes[1], marginalLikes[2]);

	
      // this can be simplified... (mf is fixed...)
      if (mf > 0.5) { // minor is unknown
	if (g == AB) {
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	} else if (g==BB) {
	  bcfgenotypes[0] = bcf_gt_unphased(true);
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	}

	if (it->genotypecall == AB) {
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	  
	  phredlikes[3] = KNOWN_PL;
	  phredlikes[4] = 0;
	  phredlikes[5] = KNOWN_PL;
	} else if (it->genotypecall == BB) {
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	  bcfgenotypes[2] = bcf_gt_unphased(true);
	  
	  phredlikes[3] = KNOWN_PL;
	  phredlikes[4] = KNOWN_PL;
	  phredlikes[5] = 0;
	} else {
	  phredlikes[3] = 0;
	  phredlikes[4] = KNOWN_PL;
	  phredlikes[5] = KNOWN_PL;
	}
	  
	logLikeToPL(marginalLikes, 3, phredlikes); // conditional PLs; the first 3

	
      } else {
	
	if (g == AB) {
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	} else if (g==BB) {
	  bcfgenotypes[2] = bcf_gt_unphased(true);
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	}
	
	if (it->genotypecall == AB) {
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	  phredlikes[0] = KNOWN_PL;
	  phredlikes[1] = 0;
	  phredlikes[2] = KNOWN_PL;
	  
	} else if (it->genotypecall == BB) {
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	  bcfgenotypes[0] = bcf_gt_unphased(true);
	    
	  phredlikes[0] = KNOWN_PL;
	  phredlikes[1] = KNOWN_PL;
	  phredlikes[2] = 0;
	} else {
	  phredlikes[0] = 0;
	  phredlikes[1] = KNOWN_PL;
	  phredlikes[2] = KNOWN_PL;
	}
	  
	logLikeToPL(marginalLikes, 3, phredlikes+3); // conditional PLs; the last 3
      }
    }



    //bcf_update_info_float(hdr, rec, "AF", &(it->af), 1);

    int32_t ad[3];
    ad[0] = counts.refQuals.size();
    ad[1] = counts.altQuals.size();
    ad[2]=bcf_int32_vector_end;
    bcf_update_info_int32(writer->hdr, rec, "AD", ad, 2);
    if (counts.otherQuals.size()) {
      int32_t dp = counts.otherQuals.size();
      bcf_update_info_int32(writer->hdr, rec, "OTH", &dp, 1);
    }
    
    
    // cast into floats
    for(i=0; i < N_GENOS; ++i)
      genolikesf[i] = genolikes[i];

    bcf_update_info_float(writer->hdr, rec, "JL", genolikesf, N_GENOS);
    
    bcf_update_genotypes(writer->hdr, rec, bcfgenotypes, 4);
    bcf_update_format_int32(writer->hdr, rec, "PL", phredlikes, 6);


    int32_t genoQuals[2];
    if (!opt.knowns || it->genotypecall==NOCALL ) {
      genoQuals[0] = getGQ(marginalLikes); // considers the three marginal likelihoods for id1
      genoQuals[1] = getGQ(&(marginalLikes[3])); // and again, for id2
    } else if (mf > 0.5) { // minor is unknown
      genoQuals[0] = getGQ(marginalLikes);
      genoQuals[1] = KNOWN_PL;
    } else {
      genoQuals[1] = getGQ(marginalLikes);
      genoQuals[0] = KNOWN_PL;
    }
	     
    bcf_update_format_int32(writer->hdr, rec, "GQ", genoQuals, 2);
    
    if (bcf_write1(writer->fp, writer->hdr, rec) ) {
      cerr << "Write error. That's not good." << endl;
      exit(1);
    }
  }

  bam_destroy1(b);
  bam_hdr_destroy(header);
  sam_close(in);
  hts_idx_destroy(bamidx);
  
  bcf_destroy1(rec);

  /*
  bcf_hdr_destroy(hdr);
  
  if (hts_close(fp))
    cerr << "Failed to close " << opt.outVCF << endl;


  if (file_format=="wb") {
    if ( bcf_index_build3(opt.outVCF.c_str(), NULL,14,opt.numThreads))
      cerr << "Some problem writing the B/VCF file index. That's not good!" << endl;
  }
  */

}


void
writeOneChrom(Options &opt, unsigned chromosome, double mf_hat, double error, VcfWriter *writer) {

  vector<Locus> loci;
  string chr = "chr" + to_string(chromosome); // autosomes only. changing this to the X will be fun....

  // grab all of the loci we want to type
  auto t1 = std::chrono::system_clock::now();
  populateLoci(&loci, opt, chr);
  
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = t2-t1;

  cerr << chr << " read : " << loci.size() << " records in " << elapsed_seconds.count() << " elapsed seconds" << endl;

  if (opt.knowns) {
    // a little clunky, but let's try it.
    // todo: check for BCF (must); not VCF.gz
    int nadded= addKnowns(loci, opt.knownBcf ,opt.knowns, chr);
    auto t3 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t3-t2;
    t2 = t3;
    cerr << chr << " addknowns " << nadded << " : elapsed seconds " << elapsed_seconds.count() << endl;
  }
  
  writeVcfOneChrom(mf_hat, error, &loci, opt, chr, writer);
  
}


// heavily influenced by: https://github.com/odelaneau/GLIMPSE/blob/master/phase/src/io/genotype_writer.cpp
// see also: http://broadinstitute.github.io/gamgee/doxygen/hts_8h.html
void
writeVcf(double mf, double error, vector<Locus> *loci, Options &opt) {


  unsigned len;
  std::string file_format = "wb";
  unsigned i;
  
  if (opt.outVCF != "-") {
    len = opt.outVCF.length();
    if (len > 4 && opt.outVCF[len-1] == 'f' &&  opt.outVCF[len-2] == 'c'  && opt.outVCF[len-3] == 'v')
      file_format = "w"; // uncompressed VCF if that's what's asked for.
  }
	      
  htsFile * fp = hts_open(opt.outVCF.c_str(),file_format.c_str());
  if (opt.numThreads > 1)
    hts_set_threads(fp, opt.numThreads);


  bcf_hdr_t *hdr = bcf_hdr_init("w");
  
  
  bcf_hdr_append(hdr, string("##source=Demixtify v" + std::string(VERSION)).c_str());
  bcf_hdr_append(hdr, string("##mixturefraction=" + std::to_string(mf)).c_str()); // yay c++11!
  bcf_hdr_append(hdr, string("##empiricalerror=" + std::to_string(error)).c_str());
  // log key items from how this file was made
  bcf_hdr_append(hdr, string("##bamFile=" + std::string(opt.bamFilename)).c_str() );
  bcf_hdr_append(hdr, string("##fstFilt=" + std::string(opt.fstTag) + ";" +  std::to_string(opt.maxFstT)).c_str() );

  // and clearly state hypotheses (along w/ reminders as to what was used)
  if (opt.knowns) {
    bcf_hdr_append(hdr, string("##hypothesis=1known+1unknown").c_str() );
    bcf_hdr_append(hdr, string("##knownIndex=" + std::to_string(opt.knowns)).c_str() );
    bcf_hdr_append(hdr, string("##knownFile=" + std::string(opt.knownBcf)).c_str() );
  } else
    bcf_hdr_append(hdr, string("##hypothesis=2unknowns").c_str() );


  
  
  // extract out the chromosomes used...
  std::string prev("");
  std::string thisc("");
  
  for (std::vector<Locus>::iterator it = loci->begin(); it != loci->end(); ++it) {
    len = it->region.length();
    
    for (i=0; i < len; ++i) {
      if (it->region[i] == ':') {
	thisc = it->region.substr(0, i);
	break;
      }
    }

    if (i==len) {
      cerr << "Should never happen " << it->region << endl;
      exit(1);
    }

    if (prev != thisc) 
      bcf_hdr_append(hdr, std::string("##contig=<ID="+ thisc + ">").c_str());
    
    prev = thisc;
  }

  bcf_hdr_append(hdr, "##FILTER=<ID=CL,Description=\"Site provides conditional likelihood\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"ALT allele frequency\">");
  bcf_hdr_append(hdr, "##INFO=<ID=JL,Number=9,Type=Float,Description=\"Joint genotype ln-likelihoods, unnormalized, in this order (Minor|Major): AA|AA,AA|AB,AA|BB,AB|AA, ... BB|BB, \">");  
  bcf_hdr_append(hdr, "##INFO=<ID=AD,Number=2,Type=Integer,Description=\"Allele Depth; Ref,Alt; filtered\">");
  bcf_hdr_append(hdr, "##INFO=<ID=OTH,Number=1,Type=Integer,Description=\"Number of non-ref/alt bases; no annotation means 0\">");
  
  bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotypes\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods; Phred scaled, marginal\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality; Phred scaled, taken as likelihood(genotype call is right)/likelihood(wrong)\">");



  bcf_hdr_add_sample(hdr, "Major");
  bcf_hdr_add_sample(hdr, "Minor");
  bcf_hdr_add_sample(hdr, NULL);      // to update internal structures ; may or may not be necessary?!
  
  if (bcf_hdr_write(fp, hdr) ) {
    cerr << "Failed to write header to file: " << opt.outVCF << endl;
    exit(1);
  }

  bcf1_t *rec = bcf_init1();

  /* data structures for converting allele counts into mixed/unmixed genotypes */
  int AAindexes[N_GENOS];
  int ABindexes[N_GENOS];
  int BBindexes[N_GENOS];

  int AAindexesMinor[N_GENOS];
  int ABindexesMinor[N_GENOS];
  int BBindexesMinor[N_GENOS];
  
  int *idx;

  double marginalLikes[N_MARGINALS];
  
  getGenoIterators("AA", MAJOR, AAindexes);
  getGenoIterators("AB", MAJOR, ABindexes);
  getGenoIterators("BB", MAJOR, BBindexes);
  
  getGenoIterators("AA", MINOR, AAindexesMinor);
  getGenoIterators("AB", MINOR, ABindexesMinor);
  getGenoIterators("BB", MINOR, BBindexesMinor);

  // convert genotypes + mixture fractions to allele weights
  double aweights[N_GENOS];
  genotypesToAlleleWeights( mf, aweights);

  double genolikes[N_GENOS];
  float genolikesf[N_GENOS];
  int phredlikes[N_MARGINALS];


  //Add records
  // 5: 2 diploid samples (4) + 1 null record (bcf_int32_vector_end)
  int bcfgenotypes[5];
  int clID = bcf_hdr_id2int(hdr,  BCF_DT_ID, "CL");

  VcfIterator knownVcf;
  if (opt.knowns) {
    knownVcf=initVcfIterator(opt.knownBcf, opt.knowns, true);
  }
  // note bedFilename is a misnomer; it's just an input filename
  VcfIterator wgsPanel = initVcfIterator(opt.bedFilename, 0, false);
  Locus loc;
  BaseCounter counts;

  // Stuff we need to process bam 
  bam1_t *b = bam_init1();
  samFile *in = sam_open(opt.bamFilename, "r");
  if (in == NULL)
    exit(EXIT_FAILURE);

  if (opt.numThreads>1)
    hts_set_threads(in, opt.numThreads);
  
  initReferenceGenome(opt.bamFilename, in, opt.referenceFasta);
      
  sam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL)
    exit(EXIT_FAILURE);
  
  hts_idx_t *bamidx;  // initiate a BAM index
  bamidx = sam_index_load(in, opt.bamFilename);   // second parameter is same as BAM file path
  if (bamidx == NULL) {
    cerr << "Failed to load index for " << opt.bamFilename << endl;
    exit(EXIT_FAILURE);
  }

  loc.af = loc.Fst=0.0;
  
  //  for (std::vector<Locus>::iterator it = loci->begin(); it != loci->end(); ++it, ++counts) {
  while ( bcf_read(wgsPanel.fp, wgsPanel.hdr, wgsPanel.rec)>=0 ) {
    bcf_unpack(wgsPanel.rec, BCF_UN_STR);

    if (!bcf_is_snp(wgsPanel.rec) || wgsPanel.rec->n_allele > 2)
      continue;
      // Clear current VCF record
    bcf_clear1(rec);
    counts.refQuals.clear();
    counts.altQuals.clear();
    counts.otherQuals.clear();
    counts.badCount=0;
    
    // todo; consider merging headers...
    string ch( bcf_hdr_id2name(wgsPanel.hdr, wgsPanel.rec->rid));    
    rec->rid = bcf_hdr_name2id(hdr, ch.c_str()); 
    rec->pos = wgsPanel.rec->pos;

    // make a locus object
    loc.region=ch + ":" + std::to_string(wgsPanel.rec->pos+1) + "-" + std::to_string(wgsPanel.rec->pos+1);
    loc.allele1 = wgsPanel.rec->d.allele[0][0];
    loc.allele2 = wgsPanel.rec->d.allele[1][0];
    loc.pos = wgsPanel.rec->pos+1;
    loc.genotypecall=NOCALL;

    // add ref and alt allele categories
    std::string alleles = std::string(1, loc.allele1) + "," + std::string(1, loc.allele2);
    bcf_update_alleles_str(hdr, rec, alleles.c_str());

    // fill in the counts data (bam processing)
    if (! summarizeRegion(in, b, header, bamidx, &loc, counts, opt)) {
      cerr << "Skipping region " << loc.region << endl; // when this happen results[i].bad == UINT_MAX
      continue;
    }

    // fill in the known genotype (if the SNP has been typed)
    if (opt.knowns) {
      loc.genotypecall=vcfGetGenotype(knownVcf, loc, opt.knowns);
    }
    
#if DEBUG
    cerr << loc <<  " " << (unsigned) loc.genotypecall << endl;
    cerr << counts << endl;
#endif
    
    // default genotypes set to 0/0 (common case)
    std::fill(bcfgenotypes, bcfgenotypes + 4, bcf_gt_unphased(false));
    bcfgenotypes[4] = bcf_int32_vector_end;
    
    // -1s turn off posterior prob calc
    // todo: use equation 3
    computeLikesWithM(&counts, aweights, error, genolikes, -1, -1, opt.trustQualities);

    if (!opt.knowns || loc.genotypecall==NOCALL ) {
      
      // corner case; minor is known, in which case the joint genotype matrix is flipped
      if (mf > 0.5)
	computeMarginals(AAindexesMinor,ABindexesMinor,BBindexesMinor,AAindexes,ABindexes,BBindexes,genolikes,marginalLikes);
      else
	computeMarginals(AAindexes,ABindexes,BBindexes,AAindexesMinor,ABindexesMinor,BBindexesMinor,genolikes,marginalLikes);
      
      //computeMarginals(AAindexes,ABindexes,BBindexes,AAindexesMinor,ABindexesMinor,BBindexesMinor,genolikes,marginalLikes);
      //logLikeToPL(const double *in, int n, int *out) {
      logLikeToPL(marginalLikes, 3, phredlikes); 
      logLikeToPL(marginalLikes+3, 3, phredlikes+3); 

      int g;
      g = LIKES_TO_GENO(marginalLikes[0], marginalLikes[1], marginalLikes[2]);
	
      if (g == AB) {
	bcfgenotypes[1] = bcf_gt_unphased(true);
      } else if (g==BB) {
	bcfgenotypes[0] = bcf_gt_unphased(true);
	bcfgenotypes[1] = bcf_gt_unphased(true);
      }

      g = LIKES_TO_GENO(marginalLikes[3], marginalLikes[4], marginalLikes[5]);      
      
      if (g == AB) {
	bcfgenotypes[3] = bcf_gt_unphased(true);
      } else if (g==BB) {
	bcfgenotypes[2] = bcf_gt_unphased(true);
	bcfgenotypes[3] = bcf_gt_unphased(true);
      }

      
    } else {

      bcf_update_filter(hdr, rec, &clID, 1);
      // and given what we know, grab the appropriate set of 3 genotypes of the joint likelihoods,
      /*
      if (mf > 0.50) {
	
	idx = AAindexes;
	if (loc.genotypecall == AB)
	  idx = ABindexes;
	else if (loc.genotypecall == BB)
	  idx = BBindexes;
	
	  
      } else {
      */
	idx = AAindexesMinor;
	if (loc.genotypecall == AB) {
	  idx = ABindexesMinor;
	} else if (loc.genotypecall == BB) {
	  idx = BBindexesMinor;
	}
	
	//}

      for (i=0; i < 3; ++i, ++idx) {
	marginalLikes[i] = genolikes[ *idx ];
      }
      
      int g = LIKES_TO_GENO(marginalLikes[0], marginalLikes[1], marginalLikes[2]);

	
      // this can be simplified... (mf is fixed...)
      if (mf > 0.5) { // minor is unknown
	if (g == AB) {
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	} else if (g==BB) {
	  bcfgenotypes[0] = bcf_gt_unphased(true);
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	}

	if (loc.genotypecall == AB) {
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	  
	  phredlikes[3] = KNOWN_PL;
	  phredlikes[4] = 0;
	  phredlikes[5] = KNOWN_PL;
	} else if (loc.genotypecall == BB) {
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	  bcfgenotypes[2] = bcf_gt_unphased(true);
	  
	  phredlikes[3] = KNOWN_PL;
	  phredlikes[4] = KNOWN_PL;
	  phredlikes[5] = 0;
	} else {
	  phredlikes[3] = 0;
	  phredlikes[4] = KNOWN_PL;
	  phredlikes[5] = KNOWN_PL;
	}
	  
	logLikeToPL(marginalLikes, 3, phredlikes); // conditional PLs; the first 3

	
      } else {
	
	if (g == AB) {
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	} else if (g==BB) {
	  bcfgenotypes[2] = bcf_gt_unphased(true);
	  bcfgenotypes[3] = bcf_gt_unphased(true);
	}
	
	if (loc.genotypecall == AB) {
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	  phredlikes[0] = KNOWN_PL;
	  phredlikes[1] = 0;
	  phredlikes[2] = KNOWN_PL;
	  
	} else if (loc.genotypecall == BB) {
	  bcfgenotypes[1] = bcf_gt_unphased(true);
	  bcfgenotypes[0] = bcf_gt_unphased(true);
	    
	  phredlikes[0] = KNOWN_PL;
	  phredlikes[1] = KNOWN_PL;
	  phredlikes[2] = 0;
	} else {
	  phredlikes[0] = 0;
	  phredlikes[1] = KNOWN_PL;
	  phredlikes[2] = KNOWN_PL;
	}
	  
	logLikeToPL(marginalLikes, 3, phredlikes+3); // conditional PLs; the last 3
      }
    }



    //bcf_update_info_float(hdr, rec, "AF", &(it->af), 1);

    int32_t ad[3];
    ad[0] = counts.refQuals.size();
    ad[1] = counts.altQuals.size();
    ad[2]=bcf_int32_vector_end;
    bcf_update_info_int32(hdr, rec, "AD", ad, 2);
    if (counts.otherQuals.size()) {
      int32_t dp = counts.otherQuals.size();
      bcf_update_info_int32(hdr, rec, "OTH", &dp, 1);
    }
    
    
    // cast into floats
    for(i=0; i < N_GENOS; ++i)
      genolikesf[i] = genolikes[i];

    bcf_update_info_float(hdr, rec, "JL", genolikesf, N_GENOS);
    
    bcf_update_genotypes(hdr, rec, bcfgenotypes, 4);
    bcf_update_format_int32(hdr, rec, "PL", phredlikes, 6);


    int32_t genoQuals[2];
    if (!opt.knowns || loc.genotypecall==NOCALL ) {
      genoQuals[0] = getGQ(marginalLikes); // considers the three marginal likelihoods for id1
      genoQuals[1] = getGQ(&(marginalLikes[3])); // and again, for id2
    } else if (mf > 0.5) { // minor is unknown
      genoQuals[0] = getGQ(marginalLikes);
      genoQuals[1] = KNOWN_PL;
    } else {
      genoQuals[1] = getGQ(marginalLikes);
      genoQuals[0] = KNOWN_PL;
    }
	     
    bcf_update_format_int32(hdr, rec, "GQ", genoQuals, 2);
    
    if (bcf_write1(fp, hdr, rec) ) {
      cerr << "Write error. That's not good." << endl;
      exit(1);
    }
  }

  bam_destroy1(b);
  bam_hdr_destroy(header);
  sam_close(in);
  hts_idx_destroy(bamidx);
  
  bcf_destroy1(rec);
  bcf_hdr_destroy(hdr);
  
  if (hts_close(fp))
    cerr << "Failed to close " << opt.outVCF << endl;


  if (file_format=="wb") {
    if ( bcf_index_build3(opt.outVCF.c_str(), NULL,14,opt.numThreads))
      cerr << "Some problem writing the B/VCF file index. That's not good!" << endl;
  }

  if (opt.knowns) {
    closeVcfIterator(knownVcf);
  }
  closeVcfIterator(wgsPanel);
}


/*

  Note: &loci is assumed to be pre-populated with the known genotype called to NOCALL (0)
 */

unsigned
addKnowns(std::vector<Locus> &loci, const char* knownBcf, int knownIndex, const std::string &chrom) {

  if (loci.size() ==0)
    return 0;
  
  VcfIterator vcf = initVcfIterator(knownBcf, knownIndex, true);

  // seek to the chromosome in question
  hts_itr_t *itr = bcf_itr_querys(vcf.idx, vcf.hdr, chrom.c_str());
  uint32_t recpos;
  char a1, a2, c1, c2;
  int ret;
  int32_t ngt, *gt_arr=NULL,ngt_arr=0;
  unsigned nloops, nrec=0, locIndex=0;
  
  // assumes diploid
  knownIndex=(knownIndex-1)*2; // genotypes stored as flattened 2D arrays 

  // .. .and this might only work for BCF (not vcf.gz). Double check
  while ((ret=bcf_itr_next(vcf.fp, itr, vcf.rec)) == 0) {
    bcf_unpack(vcf.rec, BCF_UN_STR);
    
    if (!bcf_is_snp(vcf.rec))
      continue;
    
    if (vcf.rec->n_allele > 2) {
      cerr << "Multiallelic record ignored in " << vcf.fname << endl <<
	"Please normalize (bcftools norm -m-) your VCF file!" << endl;
      continue;
    }
    
    recpos = vcf.rec->pos+1; // convert to 1-based indexing

    for ( ; locIndex < loci.size() &&recpos > loci[locIndex].pos; ++locIndex) {
      loci[locIndex].genotypecall=NOCALL;
    }

    
    if (locIndex >= loci.size())
      break;

    // 99.9% of the time, true 0 or 1 time.
    for(nloops=0 ; locIndex < loci.size() && recpos==loci[locIndex].pos; ++nloops, ++locIndex) {

      a2 = a1 =  vcf.rec->d.allele[0][0]; // by definition, must be 1 char long

      if (a1 != loci[locIndex].allele1 && a1 != loci[locIndex].allele2) {
	cerr << "Unexpected reference allele " << loci[locIndex] << " vs " << a1 << endl;
	continue;
      }
      
      // at a minimum, we have a reference basecall; we may have an alternative call too.
      if (vcf.rec->n_allele == 2) {
	a2 = vcf.rec->d.allele[1][0];	
      }

      ngt=bcf_get_genotypes(vcf.hdr, vcf.rec, &gt_arr, &ngt_arr);//FINISH
      // ensure we have genotype records (at all)
      if (ngt <= knownIndex+1) 
	continue;

      // and not missing data (haploid or diploid)
      if ( bcf_gt_is_missing( gt_arr[knownIndex]) || bcf_gt_is_missing( gt_arr[knownIndex+1]))
	continue;
      

      c1 = a1;
      if ( bcf_gt_allele(gt_arr[knownIndex]) == 1) {
	c1 = a2;
      }
      c2 = a1;
      if ( bcf_gt_allele(gt_arr[knownIndex+1]) == 1) {
	c2 = a2;
      }

      // 1/0 heterozygotes; let's convert them to 0/1
      if (c1 == loci[locIndex].allele2 && c2 == loci[locIndex].allele1) {
	c1 = loci[locIndex].allele1;
	c2 = loci[locIndex].allele2;
      }
      
      // c1,c2 are alleles (chars)
      if (c1 == loci[locIndex].allele1) { // ref matches
	
	// alt (if it exists) is consistent w/ the record
	if (c2 == loci[locIndex].allele1 || c2 == loci[locIndex].allele2) {
	  
	  if (c1 != c2) {
	    loci[locIndex].genotypecall=AB;
	  } else if (c1 == a1) { // for split multiallelic variants, bcftools norm -m - will encode 0s by default (which may get overwritten by subsequent records at the same coordinate
	    
	    if (loci[locIndex].genotypecall==NOCALL)
	      loci[locIndex].genotypecall=AA;
	    
	  }  

	  
	  ++nrec;
	}
      } else if (c1 == loci[locIndex].allele2 && c1==c2) { // alt homozygotes
	loci[locIndex].genotypecall=BB;
	++nrec;
      }
      
      
    }
    // ensures that all multiallelic records are compared against each other.
    if (nloops>1)
      locIndex -= (nloops-1);
    
    if (locIndex==loci.size())
      break;

  }

  if (gt_arr != NULL)
    free(gt_arr);
  
  closeVcfIterator(vcf);

  if (nrec==0) {
    cerr << "Failed to import any known genotype records from chromosome: " << chrom << endl << "\tReturn code is: " << ret << endl;
  }
  
#if DEBUG
  for (unsigned i=1; i < loci.size(); ++i) {
    if (loci[i-1].pos==loci[i].pos) {
      cerr << loci[i-1] << endl << loci[i] << endl;
    }
  }
  cerr << nrec << " vs " << loci.size() <<endl;

#endif
  

  if (nrec < loci.size())
    return nrec;
      
  return (unsigned) loci.size();
      
}




/*
  Augments loc , filling out the known genotype field.
 */
unsigned
addKnowns(std::vector<Locus> &loci, const char* knownBcf, int knownIndex) {


  VcfIterator v = initVcfIterator(knownBcf, knownIndex, true);

  /*
    // some simple multiallelic testing; make sure we can get either (or both) alleles
  loci[0].region="chr1:896680-896680";
  loci[0].pos=896680;
  loci[0].allele1='G';
  loci[0].allele2='T';
  loci[0].genotypecall = vcfGetGenotype(v, loci[0], knownIndex);
  cout << loci[0] << " "  << (unsigned)loci[0].genotypecall << endl;

  loci[0].pos=896680;
  loci[0].allele1='G';
  loci[0].allele2='A';
  loci[0].genotypecall = vcfGetGenotype(v, loci[0], knownIndex);
  cout << loci[0] << " "  << (unsigned)loci[0].genotypecall << endl;
  exit(1);
  */
  unsigned nrec=0;
  for (auto& loc : loci) {
    loc.genotypecall = vcfGetGenotype(v, loc, knownIndex);
    if (loc.genotypecall > 0)
      ++nrec;
    //cout << loc << " "  << (unsigned)loc.genotypecall << endl;
  }
  

  closeVcfIterator(v);
  return nrec;
}

// poor mans string compare, case insensitive.
// c can be upper or lower case, or a mix.
// c can also be a csv (, are treated as nuls;)
// b must be upper case.
// return value is the c iterated to the next word (or nul)

const char*
naiveStringCompare(const char *b, const char *c, bool *match) {
  *match=false;
  
  while (*c && *b) {
    if (! (*c==*b || *c-32==*b) )
      return c;
    ++c;
    ++b;
  }


  if ( (!*c || *c==',') && !*b) {

    *match=true;
    if (*c==',') {
      return c+1;
    }
  }
  return c;
}

char
parseMode(const char *arg1) {

  char mode=NOPE;
  
  if (arg1==NULL || *arg1==0)
    return NOPE;


  const char *c;  
  while (true) {
    bool match;

    
    c = naiveStringCompare("DESCRIBE", arg1, &match);

    if (match) {

      mode |= DESCRIBE;
      if (*c) {
	arg1=c;
	continue;
      } else
	break;

    }

    c = naiveStringCompare("DETECT", arg1, &match);
    if (match) {
      
      mode |= DETECT;
      if (*c) {
	arg1=c;
	continue;
      } else
	break;

    }
    
    c = naiveStringCompare("DEMIX", arg1, &match);
    if (match) {
      mode |= DEMIX;

      if (*c) {
	arg1=c;
	continue;
      } else
	break;
      
    }
    
    break;

  }

  
  return mode;
  
}

/*
  reads the text report (writeTextReport) created by DESCRIBE
  and extracts the mixture fraction
 */
bool
readDescribe(const char *filename, Options &opt) {

  ifstream bedFile;
  string line;
  const char *c;
  

  double curMf, curLike, bestMf, nulLike, bestLike=-std::numeric_limits<double>::max();
  opt.mixtureFraction=bestMf=curMf=-1;
  nulLike = -std::numeric_limits<double>::max();
  
  bedFile.open(filename);
  if (bedFile.fail()) 
    return false;

  
  while (getline(bedFile, line)) {

    // format is:
    //mf\t0.10\t-8003.0
    c=line.c_str();
    if (! *c)
      break;
    else if (*c != 'm' ||c[1] != 'f')
      break;

    while (*c && *c != '\t')
      ++c;

    if (! *c) {
      cerr << "Parsing failure with: " << filename << endl;
      break;
    }

    ++c; // scroll past \t
    curMf = atof(c);

    while (*c && *c != '\t')
      ++c;

    if (! *c) {
      cerr << "Parsing failure with: " << filename << endl;
      break;
    }
    
    ++c; // scroll past \t
    curLike = atof(c);

    if (bestMf==-1) { // true on the first loop thru
      bestMf=curMf;
      bestLike = nulLike = curLike;
    } else if (bestLike < curLike) {
      bestMf=curMf;
      bestLike = curLike;
    }    
  }
  bedFile.close();

  if (bestMf==-1) {
    cerr << "Parsing failure; empty file? " << filename <<endl;
    return false;
  } else if (nulLike >= bestLike) {
    cerr << "Not a mixture; use a single-source genotyping method " << filename << endl;
    return false;
  }
  
  double delta = 2*(bestLike - nulLike);
  if (delta <= CHISQ[ opt.chisqIndex ].chisq) {
    cerr << "A mixture proportion of 0 is in your confidence interval. This may not be a mixture: " << filename << endl;
  }
  
  opt.mixtureFraction=bestMf;
  return true;
}

void
writeTextReport(std::vector<MFest> &mfs, char mode, SampleStats &ss, Options &opt) {

  char buffer[1023];
  std::ostream* outfile = &std::cout;
  
  if (mode&DESCRIBE) {
    string outfileName = opt.summaryFilenameRoot + ".mf";
    if (opt.summaryFilenameRoot != "-" && opt.summaryFilenameRoot != "/dev/stdout") {

      outfile = new std::ofstream(outfileName);
    
      if (! ((std::ofstream*)outfile)->is_open()) {
	cerr << "Failed to open " << outfileName << " for writing" << endl;
	exit(1);
      }
    }
    
    for (std::vector<MFest>::iterator it = mfs.begin(); it != mfs.end(); ++it) {
      sprintf(buffer, "mf\t%0.4f\t%0.9f\n", it->mf, it->loglike); // cout's way of handling precision sucks.
      *outfile << buffer;
    }
    sprintf(buffer, "nsnps\t%f\t%d\n", ss.sequenceError, ss.nsnps);
    *outfile << buffer;

    if (outfile != &std::cout) {
      ((std::ofstream*)(outfile))->close();
      delete outfile;
    }
    
    
  }
  
  if (mode&DETECT) {

    string xtra="";
    if (mode&DESCRIBE)
      xtra="-\t"; // makes a true tsv if both modes are selected (ie, 3 columns, while the below is naturally 2 columns)
    const char *x = xtra.c_str();
    
    
    string outfileName = opt.summaryFilenameRoot + ".detect";
    if (opt.summaryFilenameRoot != "-" && opt.summaryFilenameRoot != "/dev/stdout") {
      outfile = new std::ofstream(outfileName);
    
      if (! ((std::ofstream*)outfile)->is_open()) {
	cerr << "Failed to open " << outfileName << " for writing" << endl;
	exit(1);
      }
    }
    
    sprintf(buffer, "coverage\t%s%0.4f\nnsnps\t%s%u\n", x, ss.coverage, x, ss.nsnps);
    *outfile << buffer;
    sprintf(buffer, "empirical_sequence_error\t%s%0.4f\nreported_sequencing_error\t%s%0.4f\n", x, ss.sequenceError, x, ss.sequencingError);
    *outfile << buffer;
    
    if (opt.recalibrateQualities) {
      sprintf(buffer, "initial_sequencing_error\t%s%0.4f\n", x, ss.sequencingErrorInitial);
      *outfile << buffer;
    }
    
    unsigned i, j, bestIndex=0;
    double bestAltLike, delta, bestLike=mfs[0].loglike;

    if (mfs.size() > 1) {
      bestAltLike=mfs[1].loglike;
    } else {
      bestAltLike=bestLike;
    }
    
    for ( i=1; i < mfs.size(); ++i) {
      if (mfs[i].loglike > bestLike) {
	bestLike = mfs[i].loglike;
	bestIndex=i;
      }
      if (mfs[i].loglike > bestAltLike) {
	bestAltLike=mfs[i].loglike;
      }
    }

    // for some context, see: https://en.wikipedia.org/wiki/Wilks%27_theorem
    // tl;dr, a chi square table (in some contexts) can be use to draw a confidence region around a point estimate.

    sprintf(buffer, "mp_pointestimate\t%s%0.4f\n", x, mfs[bestIndex].mf);
    *outfile << buffer;

    // log likelihood ratio
    if (mfs.size()>1) {
      double ln10=-10./PHRED_SCALAR; // equivalent to ln(10)
      sprintf(buffer, "log_10_lr\t%s%0.4f\n",  x, (bestAltLike/ln10) - (mfs[0].loglike/ln10) ); // for the forensic audience, provide the log10lr of the alt (2-person) vs null (1person) hypotheses
    }
    *outfile << buffer;
    
    i=opt.chisqIndex;
      
    // least lower bound.
    for (j=0; j < mfs.size(); ++j) {
      delta = 2*(bestLike - mfs[j].loglike);
      if (delta <= CHISQ[i].chisq) {
	break;
      }
    }
      // greatest upper bound
    sprintf(buffer, "mp_lowerbound_%0.3f\t%s%0.4f\n", CHISQ[i].percentile, x, mfs[j].mf);
    *outfile << buffer;
    
    for (j=mfs.size(); j >0; --j) {
      delta = 2*(bestLike - mfs[j-1].loglike);	
      if (delta <= CHISQ[i].chisq) {
	--j;
	break;
      }
    }
    sprintf(buffer, "mp_upperbound_%0.3f\t%s%0.4f\n", 1.0-CHISQ[i].percentile, x, mfs[j].mf);
    *outfile << buffer;

    if (outfile != &std::cout) {
      ((std::ofstream*)(outfile))->close();
      delete outfile;
    }
  }


  if (opt.recalibrateQualities) 
    return;
  
  // eg, reported base quality is 20 (average) (1% error rate), but the empirical rate is higher
  //
  if (ss.sequenceError > ss.sequencingError) {
      
    // ad hoc separation; error rate is off by more than a factor of 2
    if (ss.sequenceError > ss.sequencingError +  ss.sequencingError) {
      if (opt.trustQualities) {
	
	cerr << "Base quality warning-- the reported base qualities are systematically too large (empirical vs reported error): "
	     << ss.sequenceError << " " << ss.sequencingError  << endl;
	cerr << "Your sample may falsely present as a mixture" << endl;
      } else 
	  cerr << "Base quality warning-- the reported base qualities are very over-optimistic " << endl;
	
    
      cerr << "Strongly consider recalibrating base qualities... " << endl;
      
    } else { // bqs are off, but less than a factor of 2.
      if (opt.trustQualities) {
	cerr << "Base quality warning-- the reported base qualities appear somewhat elevated: (empirical vs reported error): "
	     << ss.sequenceError << " " << ss.sequencingError  << endl;
      } else {
	cerr << "Base quality warning-- the reported base qualities are over-optimistic; not trusting base qualities won't fix everything" << endl;
      }
    }
  }
  
  
}

int
main(int argc, char **argv) {

  Options opt;
  char mode = parseMode(argv[1]);
  if (!mode) {
    die(argv[0], "Invalid mode detected");
  }
  
  if (! parseOptions(argv, argc, opt, mode)) {
    die(argv[0], NULL);
  }
  
  
  initQual2Prob();
  
  bcf_hdr_t *header=NULL;
  vector<Locus> loci;
  BaseCounter *results=NULL;
  
  if (mode) { // should be if true.

    vector<BaseCounter> tmpResults;
    
    if (opt.parseVcf) {

      if (mode == DEMIX) {// note the == Fst tags are not necessary when we are deconvolving the mixture...
	if (NULL == (header= (bcf_hdr_t*)readPanel(opt.inBcf, loci, NULL, opt.AFtag, NULL, opt.maxFstD, opt.maxFstT)) ) {
	  die(argv[0], "Failed to parse the vcf file ... ");
	}
      } else {
	if (NULL == (header= (bcf_hdr_t*)readPanel(opt.inBcf, loci, NULL, opt.AFtag, opt.fstTag.c_str(), opt.maxFstD, opt.maxFstT)) ) {
	  die(argv[0], "Failed to parse the vcf file ... ");
	}
      }
      results = new BaseCounter[ loci.size() ];

      if (opt.knownBcf != NULL && opt.mixtureFraction < 0.) {
	unsigned nrec = addKnowns(loci, opt.knownBcf, opt.knowns);
	cerr << nrec << " known genotypes imported (from the low-fst panel)" << endl;
      }
      
    } else {

      // used for reading simulation data...
      if (! readBed(opt.bedFilename, loci, tmpResults, opt)) 
	die(argv[0], "Failed to parse the bed file ... ");

      
      if (tmpResults.size()) {

	// may happen if some records have allele counts and others do not
	if (tmpResults.size() != loci.size()) {
	  die(argv[0], "Your bed file is malformed...");
	}
	
	results=tmpResults.data(); //note: no need to free.
	opt.numThreads=0; // kludge, but we use multithreading to parse the bam file.
	// if I'm here, it means that a vector of counts was added to the bed file; hence no threads!

	if (loci[0].genotypecall != NOCALL)
	  opt.knowns=1; // treated as a flag, and an index to grab out of the vcf file
	// if we're here, there's no VCF file to parse...
      }
      
    }

    
    // borrowed from: https://dearxxj.github.io/post/5/
    // non-threaed version; open the bam file in main
    if (opt.numThreads>0) { // todo, multithreaded io? (as else)
      bam1_t *b = bam_init1();
      samFile *in = sam_open(opt.bamFilename, "r");
      if (in == NULL)
	return -1;

      initReferenceGenome(opt.bamFilename, in, opt.referenceFasta);
      
      sam_hdr_t *header = sam_hdr_read(in);
      if (header == NULL)
	return -1;
      
      hts_idx_t *idx;  // initiate a BAM index
      idx = sam_index_load(in, opt.bamFilename);   // second parameter is same as BAM file path
      if (idx == NULL) {
	cerr << "Failed to load index" << endl;
	return -1;
      }


      // this does a single seek.
      summarizeAllRegions(in, b, header, idx, loci, results, opt);

      bam_destroy1(b);
      bam_hdr_destroy(header);
      sam_close(in);
      hts_idx_destroy(idx);
      
    }

    
    SampleStats ss;
    ss.sequenceError = pow(10, opt.maxBaseQuality/-10.0);  // set to the min error (0.001)
    ss.nsnps=0;


    // bool readDescribe(const char *filename, Options &opt)
    if (opt.describeFilename!= NULL) {
      // pulls the mixture fraction from file.
      // sets opt.mixtureFraction accodingly
      if (! readDescribe(opt.describeFilename, opt) ) {
	cerr << "Parsing problem with: " << opt.describeFilename << endl <<
	  "Exiting..." << endl;
	return 1;
      }
      cout << opt.mixtureFraction << endl;
      exit(1);
    }
    
    double mf_hat = opt.mixtureFraction;

    std::vector<MFest> mfs;
    if (opt.mixtureFraction >= 0.0) // we say what the MF is. just one likelihood to assess.
      mfs.resize(1);
    else
      mfs.resize(opt.ngrid+1);

    // empirical error estiamte; truncated to be no larger than Q30 (default)
    ss.sequenceError = estimateError(results, loci.size(), ss.sequenceError);
    if (ss.sequenceError < 0) {
      cerr << "No usable data found" << endl;
      die(argv[0], NULL);
    }
    
    ss.sequencingErrorInitial=ss.sequencingError=estimateQualError(results, loci.size());
    
    // todo posthoc downsample, optionally.
    ss.coverage=estimateCoverage(results, loci.size());

        // todo add recalibration to the options.
    // also record previous quality
    if (opt.recalibrateQualities) {

      if (ss.sequenceError > ss.sequencingError) 
	recalibrateBaseQualities( (1.0- ss.sequenceError) / (1.0-ss.sequencingError) );

      ss.sequencingError=estimateQualError(results, loci.size());

    }
    
    if (opt.desiredCoverage > 0) {
      // cerr << "HERE " << ss.coverage <<  " " << 1.- opt.desiredCoverage / ss.coverage <<  endl;
      
      if (opt.desiredCoverage < ss.coverage) {
	posthocDownsample(results, loci.size(), 1. - opt.desiredCoverage / ss.coverage);
	// note we re-estimated coverage as the above is probabilistic.
	ss.coverage=estimateCoverage(results, loci.size());
	cerr << " new coverage " << ss.coverage << endl;
	
	// if the specified coverage is more than what's available...
      } else if (opt.desiredCoverage > ss.coverage + std::numeric_limits<double>::epsilon() ) {
	cerr << "Requested coverage ("<< opt.desiredCoverage <<
	  ") Exceeds the mean read depth (" << ss.coverage << ")" << endl <<
	  "No (further) downsampling will be performed" << endl;
	
      }
      
      
    }

    mf_hat = estimateMF_1Thread(results, &loci, &opt, ss.sequenceError, mfs, &ss.nsnps);

    if (mode&DETECT || mode&DESCRIBE) {
      writeTextReport(mfs, mode, ss, opt);
    }

    if (! (mode&DEMIX) )
      exit(0);

    
      
    VcfWriter writer;
    writer.fp=NULL;
    writer.hdr=NULL;
    writer.results=results; // these are incremented when writing (todo)
    writer.loci = &loci;
    
    for (unsigned i =1 ; i < 23; ++i) {
      // todo; call writeVcfOneChrom instead; control results/loci accordingly. Note that writer has both.
      writeOneChrom(opt, i, mf_hat, ss.sequenceError, &writer);
    }
    
    closeVcfWriter(writer, opt);

    
    if (!tmpResults.size())
      delete[] results;
    
  }

  // if writing a BCF/VCF, I need to recycle the header
  if (header != NULL) {
    bcf_hdr_destroy(header);    
  }
  
  return EXIT_SUCCESS;
}


// takes in a filename,
// initializes a VcfIterator
// exits on failure.
VcfIterator
initVcfIterator(const char *fname, int knownIndex, bool needidx) {
  VcfIterator vcf;
  vcf.fp    = hts_open(fname,"r");
  if (vcf.fp==NULL) {
    cerr << "Failed to open " << fname << endl;
    exit(1);
  }

  vcf.hdr = bcf_hdr_read(vcf.fp);
  
  if (vcf.hdr==NULL) {
    cerr << "Failed to parse the header of " << fname << endl;
    exit(EXIT_FAILURE);
  }


  if (needidx) {
    vcf.idx = bcf_index_load(fname);
    if (vcf.idx == NULL) {
      cerr << "Failed to open the index to: " << fname << endl;
      exit(EXIT_FAILURE);
    }
  } else
    vcf.idx=NULL;
  
  vcf.rec    = bcf_init();
  vcf.fname=fname;

  if (knownIndex > 0) {
    if (knownIndex > bcf_hdr_nsamples(vcf.hdr)) {
      cerr << endl << "You asked for sample number: " << knownIndex << ", but your VCF file doesn't have that many samples!" << endl << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  return vcf;
}

bool
closeVcfIterator(VcfIterator &vcf) {

  bcf_destroy(vcf.rec);
  bcf_hdr_destroy(vcf.hdr);


  if (vcf.idx != NULL)
    hts_idx_destroy(vcf.idx);

  if ( hts_close(vcf.fp) ) {
    cerr << "hts_close(" << vcf.fname << ") non-zero status" << endl;
    return false;
  }
  
  return true;
}

/*
Takes in a locus (loc)
a vcf record (vcf), 
and the index of the known contributor (typically 0)
and returns the true genotype as:
#define NOCALL 0
#define AA 1
#define AB 2
#define BB 3
 */
char
vcfGetGenotype(VcfIterator &vcf, Locus &loc, int knownIndex) {
  hts_itr_t *itr = bcf_itr_querys(vcf.idx, vcf.hdr, loc.region.c_str());
  char genotype = NOCALL;
  uint32_t recpos;
  char a1, a2;
  int ret;
  int32_t ngt, *gt_arr=NULL,ngt_arr=0;

  // assumes diploid
  knownIndex=(knownIndex-1)*2; // genotypes stored as flattened 2D arrays 

  // .. .and this might only work for BCF (not vcf.gz). Double check
  while ((ret=bcf_itr_next(vcf.fp, itr, vcf.rec)) == 0) {
    bcf_unpack(vcf.rec, BCF_UN_STR);

    
    if (!bcf_is_snp(vcf.rec))
      continue;
    
    if (vcf.rec->n_allele > 2) {
      cerr << "Multiallelic record ignored in " << vcf.fname << endl <<
	"Please normalize (bcftools norm -m-) your VCF file!" << endl;
      continue;
    }
    
    recpos = vcf.rec->pos+1; // convert to 1-based indexing

    //cerr << bcf_hdr_id2name(vcf.hdr, vcf.rec->rid) <<  " " << loc.pos << " " << recpos << endl;
    // we looked, but got too far!
    if (loc.pos < recpos)
      break;

    int chrpos = loc.region.find_first_of(':');
    if (chrpos < 0) {
      cerr << "Parse error with locus " << loc << endl;
      exit(EXIT_FAILURE);
    }
    
    // todo: check chromosome...
    int chrid =  bcf_hdr_name2id(vcf.hdr, loc.region.substr(0, chrpos).c_str());
    if (chrid !=  vcf.rec->rid ) {
      // incredibly unlikely, but we may scroll past the end of the target chromosome?
      break;
    }

    if (loc.pos == recpos) {
      a1 =  vcf.rec->d.allele[0][0]; // by definition, must be 1 char long

      if (a1 != loc.allele1 && a1 != loc.allele2) {
	cerr << "Unexpected reference allele " << loc << " vs " << a1 << endl;
	continue;
      }
      
      // at a minimum, we have a reference basecall; we may have an alternative call too.
      if (vcf.rec->n_allele == 2) {
	a2 = vcf.rec->d.allele[1][0];
	
	// can happen (multiallelic site; normalized)
	if (a2 != loc.allele1 && a2 != loc.allele2)
	  continue;
	
      }

      ngt=bcf_get_genotypes(vcf.hdr, vcf.rec, &gt_arr, &ngt_arr);
      if (ngt <= knownIndex+1) 
	break;

      if ( bcf_gt_is_missing( gt_arr[knownIndex]) || bcf_gt_is_missing( gt_arr[knownIndex+1]))
	break;

      
      a1 = bcf_gt_allele(gt_arr[knownIndex]);
      a2 = bcf_gt_allele(gt_arr[knownIndex+1]);
      if (a1 != a2) 
	genotype=AB;
      else if (a1 == 0)
	genotype=AA;
      else if (a1 == 1)
	genotype=BB;
      else // this should not happen. (at all. ever)
	cerr << "Unexpected genotype " << (unsigned) a1 << " " << (unsigned) a2 << " " << endl << loc << endl;
	
      

      break;
    }
    //printf("ret: %d - id: %d, pos: %d, len: %d\n", ret, rec->rid, rec->pos, rec->rlen);

  }
  //cerr << "Return " << ret << endl;
  
  bcf_itr_destroy(itr);
  if (gt_arr != NULL)
    free(gt_arr);
  
  return genotype;
}

// This is used to estimate the MF. To do so, 
// the entirety of this VCF is read into memory.
// adapted from: readVcf.c : https://gist.github.com/sujunhao
// input files may either be of type bed (uncompressed) or
// something that hts_open can read (bcf/vcf/vcf.gz would make sense)
void *
readPanel(const char *fname, vector<Locus> &loci, BaseCounter *results, const char* AFtag, const char* fstTag, double maxFstD, double maxFstT) {

  htsFile *fp    = hts_open(fname,"r");
  if (fp==NULL) {
    cerr << "Failed to open " << fname << endl;
    return NULL;
  }
  
  //read header
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec    = bcf_init();

  string chrom;
  int curRid=-1;

  float* af=NULL;
  int naf=0;
  
  //save for each vcf record
  while ( bcf_read(fp, hdr, rec)>=0 ) {
      //unpack for read REF,ALT,INFO,etc 
    bcf_unpack(rec, BCF_UN_STR);
    bcf_unpack(rec, BCF_UN_INFO);

    if (bcf_is_snp(rec) && rec->n_allele == 2) {
      Locus loc;      
      // returns >= on success..
      if (bcf_get_info_float(hdr, rec, AFtag, &af, &naf) < 0) {
	continue;
      } else {
	loc.af= af[0];
      }


      if (fstTag==NULL || maxFstT <= 0.) {
	loc.Fst=0.;
	// and now grab the Fst tag (should it exist)
      } else if (bcf_get_info_float(hdr, rec, fstTag, &af, &naf) < 0) {
	continue;
      } else {
	loc.Fst= af[0];
	if (loc.Fst >  maxFstD)
	  continue;
	else if (loc.Fst > maxFstT)
	  loc.Fst=maxFstT;

      }
      
      // convert the internal (int-based) chromosome representation into its string equivalent
      if (rec->rid != curRid) {
	string ch( bcf_hdr_id2name(hdr, rec->rid) );
	chrom = ch;
	curRid = rec->rid;
	
	bool gotNum=false;
	for (auto &c : ch) {
	  if (c >= '0' && c <= '9') {
	    gotNum=true;
	    break;
	  }
	}
	// forcibly only consider the autosomes
	// any sane sorting of the VCF file should put the numeric chromosomes (ie, the canonical autosomes) first
	if (! gotNum)
	  break;
      }

      //rec->pos is 0-based; convert to 1-based
      loc.region = chrom + ':' +  to_string(rec->pos+1) + '-' +  to_string(rec->pos+1);
      loc.pos = rec->pos+1;
      loc.allele1 = rec->d.allele[0][0]; // by definition, must be 1 char long
      loc.allele2 = rec->d.allele[1][0];
      loc.genotypecall=NOCALL;

      
      loci.push_back(loc);
    }

  }
    
  bcf_destroy(rec);
  // bcf_hdr_destroy(hdr); // must be freed later. TODO!

  if (af!=NULL)
    free(af);
  
  int ret;
  
  if ( (ret=hts_close(fp)) ) {
    cerr << "hts_close(" << fname << "): non-zero status" << endl;
    return NULL;
  }
  return hdr;
}

bool
readBed(char *filename, vector<Locus> &loci, vector<BaseCounter> &bc, const Options &opt) {
  ifstream bedFile;
  string line;
  
  bedFile.open(filename);
  if (bedFile.fail()) 
    return false;

  bool ret=true;
  bool printWarn=false;
  std::vector<double> afs;
  int nrecs;

  unsigned refCount, altCount, i;
  
  while (getline(bedFile, line)) {
    if (line[0] == '#') // bed files may have a header; headers must have a #. technically these can be anywhere in the file, but, well, laziness :)
      continue;
    
    stringstream ss(line);
    vector<std::string> records;
    string tmp;
    while(getline(ss, tmp, '\t')) {
      records.push_back(tmp);
    }

    if (!records.size()) // blank lines are A okay
      continue;
    
    if (records.size() < 5) {
      cerr << "Error on line: " << endl << line << endl;
      ret=false;
      break;
    }
		
    int startPosition = atoi(records[1].c_str()) + 1;
    int stopPosition=atoi(records[2].c_str());
    if (startPosition<=1 || startPosition>stopPosition) {
      cerr << "Invalid positions in line: " << line << endl;
      ret = false;
      break;
    }
    
    if (records[3].length() != 1 || records[4].length() != 1) {
      cerr << "Alleles must be of length exactly one..." << endl;
      ret = false;
      break;
    }
    
    Locus loc;
    loc.region = records[0] + ":" + to_string(startPosition) + "-" + records[2];
    loc.pos = startPosition;
    loc.allele1=records[3][0];
    loc.allele2=records[4][0];
    loc.genotypecall=NOCALL;
    loc.af=DEFAULT_AF;
    loc.Fst=0.;
    
    // convert to upper-case letters.
    if (loc.allele1 > 'Z') {
      loc.allele1 -= 32;
    }
    if (loc.allele2 > 'Z') {
      loc.allele2 -= 32;
    }
    if (!validNuc(loc.allele1) || !validNuc(loc.allele2) ) {
      if (!printWarn) {
	cerr << "Allelic states may only be A,C,G or T" << endl;
	printWarn=true;
      }
      continue;
    } else if (loc.allele1==loc.allele2) {
      cerr << "Alleles may not be the same... " << endl << line << endl;
      continue;
    }

    // it's inelegant, but for testing purposes it can be nice to include raw counts
    // note that data written this way cannot be exported to the VCF file format (eg, the genome version is not known)
    if (records.size() > 7) {
      BaseCounter foo;
      // changed in version 0.2
      // we now record the quality scores associated w/ the reference and alternative alleles.
      refCount=atoi(records[5].c_str());
      foo.refQuals.reserve(refCount);
      for (i=0; i < refCount; ++i)
	foo.refQuals.push_back( opt.maxBaseQuality);
      
      altCount=atoi(records[6].c_str());
      foo.altQuals.reserve(altCount);
      for (i=0; i < altCount; ++i)
	foo.altQuals.push_back( opt.maxBaseQuality);

      
      unsigned otherCount=atoi(records[7].c_str());
      foo.otherQuals.reserve(otherCount);
      for (i=0; i < otherCount; ++i)
	foo.otherQuals.push_back( opt.maxBaseQuality);
      
      foo.badCount=0;
      bc.push_back(foo);
    } 
    
    if (records.size() > 8) { // modified to have the allele frequency column be a csv; the 0th value is "the" allele frequency
      // any subsequent values are population-specific allele frequencies
      afs.clear();
      nrecs=csv2vec(records[8].c_str(), &afs);
      if (nrecs<1) {
	cerr << "Should never happen. " << records[8] << endl;
	exit(1);
      }
      //loc.af = atof(records[8].c_str());
      loc.af = afs[0];
      
      if (loc.af < 0. || loc.af > 1.) {
	cerr << "Invalid allele frequency on line" << endl << line << endl << records[7] << endl;
	return false;
      }

      // default Fst is 0; equivalent to no theta correction
      if (afs.size() > 1) {

	loc.Fst = afs[1];
	if (loc.Fst > opt.maxFstD) 
	  bc.back().badCount=IGNORE_SNP;
	else if (loc.Fst > opt.maxFstT)
	  loc.Fst = opt.maxFstT;
	
      }

    }
    // assumes a 0/1/2 encoding for one of the "known" genotypes
    if (records.size() > 9) {
      // we encode using two bits, 0 meaning I don't know, then 1, 2, 3 for AA, AB, BB
      // the simulation program uses 1.5 bits (0,1,2)==(AA,AB,BB). hence the +1
      loc.genotypecall = 1 + atoi(records[9].c_str());
      if (loc.genotypecall < 0 || loc.genotypecall > BB) {
	cerr << "Invalid genotype on line" << endl << line << endl << records[8] << endl;
	return false;
      }

    }
    
    loci.push_back(loc);
  }

  
  bedFile.close();
  return ret;
}


/*
  Used to write the raw counts to file.
  This is run iff opt.outCounts is not NULL
 */
void
writeCounts(const char *outfile, vector<Locus> *loci, BaseCounter *counts, double mf, double error, const Options *opt) {

  ofstream fh;
  fh.open(outfile);
  int i;

  int AAindexes[N_GENOS];
  int ABindexes[N_GENOS];
  int BBindexes[N_GENOS];

  int AAindexesMinor[N_GENOS];
  int ABindexesMinor[N_GENOS];
  int BBindexesMinor[N_GENOS];
  
  int *idx;

  double marginalLikes[N_MARGINALS];
  
  
    // nIndexes==3 in all cases.
    // these provide the integer offsets into the 9 genotype likelihoods that are associated w/ each proposed (marginal genotype), beginning with the major being AA
  getGenoIterators("AA", MAJOR, AAindexes);
  getGenoIterators("AB", MAJOR, ABindexes);
  getGenoIterators("BB", MAJOR, BBindexes);
  
  getGenoIterators("AA", MINOR, AAindexesMinor);
  getGenoIterators("AB", MINOR, ABindexesMinor);
  getGenoIterators("BB", MINOR, BBindexesMinor);


  
  // convert genotypes + mixture fractions to allele weights
  double aweights[N_GENOS];
  genotypesToAlleleWeights( mf, aweights);

  double genolikes[N_GENOS];
  //  int phredlikes[N_MARGINALS];
  
  if (fh.is_open()) {
    fh << "Region\tA\tB\tMixFrac\tErrorHat\tRefCount\tAltCount\tOtherCount\tBadCount";

    for (i =0; i < N_GENOS; ++i) {
      fh << "\t" << POSSIBLE_GENOS[i];
    }

    fh << "\tAA\tAB\tBB\tAAminor\tABminor\tBBminor";
    
    fh << endl;
    
    for (std::vector<Locus>::iterator it = loci->begin(); it != loci->end(); ++it, ++counts) {

      if (counts->badCount ==SKIP_SNP)
	continue;
      
      fh << it->region << "\t" << it->allele1 << "\t" << it->allele2 << "\t" << mf << "\t" << error <<
	"\t"  << counts->refQuals.size() <<
	"\t" << counts->altQuals.size() <<
	"\t" << counts->otherQuals.size() <<
	"\t" << counts->badCount;

      // note: -1 forces us to just compute the likelihood... (and not the posterior)
      computeLikesWithM(counts, aweights, error, genolikes, -1, -1, opt->trustQualities);
      //logLikeToPL(genolikes, N_GENOS, phredlikes);
      for (i =0; i < N_GENOS; ++i) {
	fh << "\t" << genolikes[i];
	//fh << "\t" << (int) phredlikes[i];
      }
      if (!opt->knowns || it->genotypecall==NOCALL ) {

	// okay, this is wrong. 
	if (mf > 0.5)
	  computeMarginals(AAindexesMinor,ABindexesMinor,BBindexesMinor,AAindexes,ABindexes,BBindexes,genolikes,marginalLikes);
	else
	  computeMarginals(AAindexes,ABindexes,BBindexes,AAindexesMinor,ABindexesMinor,BBindexesMinor,genolikes,marginalLikes);
	
	for (i =0; i < N_MARGINALS; ++i)
	  fh << "\t" << marginalLikes[i];
	
      } else { // todo: case where "known" is a missing genotype call; make it tabular as per the 2 unknown case

	// and given what we know, grab the appropriate set of 3 genotypes of the joint likelihoods,

	if (mf > 0.50) {
	
	  idx = AAindexes;
	  if (it->genotypecall == AB)
	    idx = ABindexes;
	  else if (it->genotypecall == BB)
	    idx = BBindexes;

	// and print them!
	  for (i=0; i < 3; ++i, ++idx) {
	    fh << "\t" << genolikes[ *idx ];
	  }
	  
	  if (it->genotypecall == AB) {
	    fh << "\t-10\t0\t-10";
	  } else if (it->genotypecall == BB) {
	    fh << "\t-10\t-10\t0";
	  } else {
	    fh << "\t0\t-10\t-10";
	  }
	  
	} else {

	  // hot encode the known genotype...
	  idx = AAindexesMinor;
	  if (it->genotypecall == AB) {
	    idx = ABindexesMinor;
	    fh << "\t-10\t0\t-10";
	  } else if (it->genotypecall == BB) {
	    idx = BBindexesMinor;
	    fh << "\t-10\t-10\t0";
	  } else {
	    fh << "\t0\t-10\t-10";
	  }
	// and print them!
	  for (i=0; i < 3; ++i, ++idx) {
	    fh << "\t" << genolikes[ *idx ];
	  }

	}
	
      }
      
      fh << endl;
    }
    fh.close();
  } else {
    cerr << "Failed to open the counts file" << outfile << endl;
  }

}
