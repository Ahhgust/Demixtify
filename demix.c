#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <threads.h>

// To future travelers:
// this presumes that there is a local version of htslib in the PWD
// feel free to modify the headers below as needed
#include "htslib/htslib/sam.h"
#include "htslib/htslib/hts.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/kseq.h"

#include "demix.h"

// sensible defaults (I hope!)
#define DEFAULT_MIN_MAPQ 20
#define DEFAULT_MIN_BASEQ 20
#define DEFAULT_MIN_READLEN 30
#define DEFAULT_NGRID 100

// used to exclude reads on the basis of the following:
#define DEFAULT_READ_FILTER (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)

// additionally, reads also must have ALL the following bitfields set
// (ie, must not trip the above filters, and must also be a proper pair)
#define DEFAULT_READ_INCLUDE_FILTER (BAM_FPROPER_PAIR)

using namespace std;


char BUFF[1024];

// make it so I can print a locus...
ostream& operator<<(ostream &os, const Locus &loc) {
  return os << loc.region << "\t" << loc.allele1 << "/" << loc.allele2;
}


void
die(const char *arg0, const char *extra) {
	
  if (extra != NULL) {
    cerr << extra << endl;
  }
  
  cerr << "Usage " << endl <<
    arg0 << " -b bamFile -d bedFile OR" << endl <<
    arg0 << " -b bamFile -v (b/v)cfFile OR" << endl <<
    arg0 << " -b bamFile -r UCSC-style-region -1 allele1 -2 allele2" << endl <<
    
    "Options " << endl <<

    "\t-h (prints this message) " << endl <<
    "\t-c countsFile (writes the allele counts to file)" << endl <<
    "\t-o outFile (writes deconvolved data in the V/BCF file format; defaults to standard output)" << endl <<
    "\t-i (includes basecalls that are adjacent to an indel). Defaults to excluding such sites" << endl <<
    "\t-g grid_size (for the grid-search on mixture fraction; the size of the grid (i.e., the precision)). Default: " << DEFAULT_NGRID  << endl <<
    "\t-k knownIndex (from the VCF file, the 1-based index of the known contributor; defaults to 0 (no known contributor))" << endl<<
    "\t-t nthreads (number of threads to use. Defaults to 1) " << endl <<    
    "\t-m mapping_quality (minimum mapping quality). Default: " << DEFAULT_MIN_MAPQ <<  endl <<
    "\t-F fraction (the mixture fraction; if unspecified, this is estimated)" << endl <<
    "\t-q base_quality (minimum base quality). Default: " << DEFAULT_MIN_BASEQ << endl <<
    "\t-L length (minimum read length). Default: " << DEFAULT_MIN_READLEN << endl <<
    "\t-f read_filter (excludes reads according to SAM read filter flags). Default: 0x" << std::hex << DEFAULT_READ_FILTER << endl <<
    "\t-g read_include_filter (includes reads if all filters are met). Default: 0x" << std::hex << DEFAULT_READ_INCLUDE_FILTER  << endl << 
    "\t-r region (genomic region, UCSC-style)" << endl <<
    "\t-1 allele1 (for -r, the first allele)" << endl <<
    "\t-2 allele2 (for -r, the second allele)" << endl << endl <<
    "For either of the filter options (-f/-g) see:" << endl <<
    "https://broadinstitute.github.io/picard/explain-flags.html" << endl << endl;

  exit(EXIT_FAILURE);
}

bool
validNuc(char c) {
  if (c == 'A' || c== 'C' || c == 'G' || c == 'T')
    return true;
  return false;
}

bool
parseOptions(char **argv, int argc, Options &opt, Locus &loc) {
	int i = 1;
	opt.filter=DEFAULT_READ_FILTER;//(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
	opt.include_filter=DEFAULT_READ_INCLUDE_FILTER;
	opt.outCounts=opt.bedFilename = opt.bamFilename= NULL;
	opt.minBaseQuality= DEFAULT_MIN_BASEQ;
	opt.minMapQuality=DEFAULT_MIN_MAPQ;
	opt.minReadLength=DEFAULT_MIN_READLEN;
	opt.outVCF="-";
	opt.knowns=0;
	
	opt.filterIndelAdjacent=true;
	opt.help=false;
	opt.parseVcf=false;
	opt.ngrid=DEFAULT_NGRID;
	opt.numThreads=1;
	opt.mixtureFraction=-1.0;
	
	loc.allele1 = loc.allele2 = 'A';
	loc.region="";
	
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
	  } else { // option has an argument
	    ++i;
	    if (i==argc) {
	      cerr << "Flag " << flag << " expects and argument!" << endl;
	      ++errors;
	      break;
	    }	
	    if (flag == 'q') {
	      opt.minBaseQuality = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'm') {
	      opt.minMapQuality = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'k') {
	      opt.knowns = max(atoi(argv[i]), 0); // make this unsigned-friendly
	    } else if (flag == 'L') {
	      opt.minReadLength = max(atoi(argv[i]), 0);
	    } else if (flag == 'f') {
	      opt.filter = atoi(argv[i]);
	    } else if (flag == 'g') {
	      opt.include_filter = atoi(argv[i]);
	    } else if (flag == 'F') {
	      opt.mixtureFraction = max(atof(argv[i]), 0.0000001);
	    } else if (flag == 't') {
	      opt.numThreads = max(atoi(argv[i]), 1);
	    } else if (flag == 'g') {
	      opt.ngrid = max(atoi(argv[i]), 2);
	    } else if (flag == 'r') {
	      loc.region = string(argv[i]);
	    } else if (flag == 'o') {
	      opt.outVCF = string(argv[i]);
	    } else if (flag == 'b') {
	      opt.bamFilename = argv[i];
	    } else if (flag == 'g') {
	      opt.ngrid = max(atoi(argv[i]), 0);
	    } else if (flag == 'd') {
	      opt.bedFilename = argv[i];
	    } else if (flag == 'c') {
	      opt.outCounts = argv[i];
	    } else if (flag == 'v') {
	      opt.bedFilename = argv[i];
	      opt.parseVcf=true;
	    } else if (flag == '1') {
	      loc.allele1 = argv[i][0];
	      if (argv[i][1] != 0) {
		cerr << "Only single-nucleotide regions are allowed" << endl;
		++errors;
	      }
	    } else if (flag == '2') {
	      loc.allele2 = argv[i][0];
	      if (argv[i][1] != 0) {
		cerr << "Only single-nucleotide regions are allowed" << endl;
		++errors;
	      }
	      
	    } else {
	      ++errors;
	      cerr << "Unexpected flag: " << argv[i-1] << endl;
	    }
	  }
		
	}
	
	if (opt.bamFilename==NULL) {
	  ++errors;
	  cerr << "No bam file was specified... I kinda need one of those!" << endl;
	}
	
	if (opt.bedFilename==NULL) {
	  if (!loc.region.length()) {
	    ++errors;
	    cerr << "No bed file was specified and no region was specified... I kinda need one of those!" << endl;		
	  }
	  // convert to upper-case letters.
	  if (loc.allele1 > 'Z') {
	    loc.allele1 -= 32;
	  }
	  if (loc.allele2 > 'Z') {
	    loc.allele2 -= 32;
	  }
	  
	  if (!validNuc(loc.allele1) || !validNuc(loc.allele2) ) {
	    ++errors;
	    cerr << "Allelic states may only be A,C,G or T" << endl;
	  } else if (loc.allele1==loc.allele2) {
	    ++errors;
	    cerr << "And alleles may not be the same... " << endl;
	  }
	}
	
	
	
	if (errors)
	  return false;


	if (opt.help)
	  die(argv[0], NULL);
	
	return true;
}


ostream& operator<<(ostream &os, const BaseCounter &b) {
  return os << "Ref\t" << b.refCount << "\tAlt\t" << b.altCount << "\tOther\t" << b.otherCount << "\tBad\t" << b.badCount;
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
    buf[i]="=ACMGRSVTWYHKDBN"[bam_seqi(s,
				       i)];
  }
  buf[i] = 0;
}


bool
summarizeRegion(samFile *in, bam1_t *b, sam_hdr_t *header, hts_idx_t *idx, Locus *loc, BaseCounter &count, Options &opt) {
  
  int coord, unused;
  hts_itr_t *iter;    // initiate an iterator
  iter  = sam_itr_querys(idx, header, loc->region.c_str());
  
  memset(&count, 0, sizeof(BaseCounter));

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
  
  while ( sam_itr_next(in, iter, b) >= 0) {

    bam1_core_t core = b->core;
    // for paired-end/mate paired. assume mates/pairs have the same name
    string readname( bam_get_qname(b) );


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
	  basecall.q = qual[offset+readPos];
	  
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
    }
	
    ++tot;
  }
// evaluate the pileup (correctly accounting for overlapping reads)
  for (auto iter = readPile.begin(); iter != readPile.end(); ++iter) {
    BQ basecall = iter->second;
    if (basecall.b==loc->allele1)
      ++count.refCount;
    else if (basecall.b==loc->allele2)
      ++count.altCount;
    else
      ++count.otherCount;
  }

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
inline void
computeLikesWithM(int acount, int bcount, double *w, double e, double out[N_GENOS]) {


  int i = 0;
  for ( ; i < N_GENOS; ++i, ++out, ++w) {
    *out = LOG_LIKE(acount, bcount, *w, e);
  }

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

double
estimateMF_1Thread(const BaseCounter *counts,  vector<Locus> *loci, const Options *opt) {

  int nLoci = (int) loci->size();
  int i,j, gridSize = opt->ngrid;
  double mf; // the proposed mixture fraction
  double error= pow(10, opt->minBaseQuality/-10.0);


  double loglike;
  double bestLike=-1;
  double bestMF=-1;
  double gridSizeD=gridSize;
  
  double aweights[N_GENOS];
  //auto genolikes = new double[ nLoci ][N_GENOS];
  double genolikes[N_GENOS];
  const BaseCounter *c;
  std::vector<Locus>::iterator it; 

  int AAindexes[N_GENOS];
  int ABindexes[N_GENOS];
  int BBindexes[N_GENOS];
  int *idx;
  int nIndexes=0;
  
  bool haveKnowns = opt->knowns>0;

  if (haveKnowns) {
    // nIndexes==5 in all cases.
    nIndexes = getGenoIterators("AA", EITHER, AAindexes);
    getGenoIterators("AB", EITHER, ABindexes);
    getGenoIterators("BB", EITHER, BBindexes);

  } else { // both unknowns, then the mixture fraction only meaningfully varies from 0-0.5 (or from 0.5-1, take you pick!)
    gridSize /=2;
    gridSizeD = gridSize;
  }
  
  // uniformly select points in the range [0,1] for the mixture fraction
  for (i = 0; i <= gridSize; ++i) {
    mf = i/gridSizeD;
    c = counts;

    loglike=0.;

    // convert the 9 genotype calls with the proposed mixture fraction
    // to get the fraction of reads expected of an "A" allele (a continuous "genotype")
    genotypesToAlleleWeights( mf, aweights);
    it = loci->begin(); 
    for (j=0; j < nLoci; ++j, ++c, ++it) {
      
      if (c->badCount !=SKIP_SNP) {
	computeLikesWithM(c->refCount, c->altCount, aweights, error, genolikes);
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
	  
	} else 
	  loglike += log_sum_exp(genolikes);
      }
    }

    cout << mf << "\t" << loglike << endl;
    
    if (i==0 || bestLike < loglike) {
      bestLike=loglike;
      bestMF=mf;
    }

  }
  
  
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

// TODO!
void
writeVcf(double mf, BaseCounter *results, vector<Locus> *loci, bcf_hdr_t *hdr) {


}

int
main(int argc, char **argv) {


  Options opt;
  Locus loc;
  if (! parseOptions(argv, argc, opt, loc)) {
    die(argv[0], NULL);
  }

  bcf_hdr_t *header=NULL;
  

  // used for testing purposes almost exclusively; just prints the base-count stats of 1 locus.
  if (loc.region.length()) {
    BaseCounter result;
    // borrowed from: https://dearxxj.github.io/post/5/
    // non-threaed version; open the bam file in main
    bam1_t *b = bam_init1();
    samFile *in = sam_open(opt.bamFilename, "r");
    if (in == NULL)
      return -1;
    
    sam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL)
      return -1;
    
    hts_idx_t *idx;  // initiate a BAM index
    idx = sam_index_load(in, opt.bamFilename);   // second parameter is same as BAM file path
    if (idx == NULL) {
      cerr << "Failed to load index" << endl;
      return -1;
    }
    if (! summarizeRegion(in, b, header, idx, &loc, result, opt)) {
      cerr << "Skipping region " << loc.region << endl;
      return -1;
    }
    cout << result << endl;
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    hts_idx_destroy(idx);
    exit(EXIT_SUCCESS);
    
  } else {
    vector<Locus> loci;

    if (opt.parseVcf) {
      if (NULL == (header= (bcf_hdr_t*)readVcf(opt.bedFilename, loci, opt.knowns)) ) {
	die(argv[0], "Failed to parse the vcf file ... ");
      }
    } else if (! readBed(opt.bedFilename, loci)) {
      die(argv[0], "Failed to parse the bed file ... ");
    }

    BaseCounter *results = new BaseCounter[ loci.size() ];

    if (opt.numThreads==1) {
      int i=0;
      // borrowed from: https://dearxxj.github.io/post/5/
      // non-threaed version; open the bam file in main
      bam1_t *b = bam_init1();
      samFile *in = sam_open(opt.bamFilename, "r");
      if (in == NULL)
	return -1;
      
      sam_hdr_t *header = sam_hdr_read(in);
      if (header == NULL)
	return -1;
      
      hts_idx_t *idx;  // initiate a BAM index
      idx = sam_index_load(in, opt.bamFilename);   // second parameter is same as BAM file path
      if (idx == NULL) {
	cerr << "Failed to load index" << endl;
	return -1;
      }
  

      
      for (std::vector<Locus>::iterator it = loci.begin(); it != loci.end(); ++it, ++i) {

	if (! summarizeRegion(in, b, header, idx, &(*it), results[i], opt)) {
	  cerr << "Skipping region " << loc.region << endl; // when this happen results[i].bad == UINT_MAX
	}
      }

      bam_destroy1(b);
      bam_hdr_destroy(header);
      sam_close(in);
      hts_idx_destroy(idx);


      
    } else {
      // multithreaded version; each
      // thread has its own file handle
      // (bam IO has state)
      
      // helps later on; make sure the number of threads <= the number of loci
      opt.numThreads = min(opt.numThreads, (int)loci.size());
	    
      thrd_t *threads = new thrd_t[ opt.numThreads ];

      SummarizeRegionHelper  *helpers = new SummarizeRegionHelper[ opt.numThreads ];
      for (int i = 0; i < opt.numThreads; ++i) {
	helpers[i] = {&loci, results, opt, i};
	thrd_create(&threads[i], threadSummarizeRegion, &helpers[i]);
      }

      
      for (int i = 0; i < opt.numThreads; ++i) {
	thrd_join(threads[i], NULL);
      }
      
      delete[] threads;
      delete[] helpers;
      // TODO (?) multithreaded version. (maybe).

    }

    double mf_hat = estimateMF_1Thread(results, &loci, &opt);
    //    deconvolveSample(    
    if (opt.outCounts != NULL) {
      writeCounts(opt.outCounts, &loci, results, mf_hat, &opt);
    }

    delete[] results;
  }

  // if writing a BCF/VCF, I need to recycle the header
  if (header != NULL) {
    bcf_hdr_destroy(header);    
  }
  
  return EXIT_SUCCESS;
}

// adapted from: readVcf.c : https://gist.github.com/sujunhao
// input files may either be of type bed (uncompressed) or
// something that hts_open can read (bcf/vcf/vcf.gz would make sense)
void *
readVcf(char *fname, vector<Locus> &loci, int knownIndex) {
	//open vcf file
	// todo:check bcf, vcf.gz, wrt to rb
  htsFile *fp    = hts_open(fname,"r");
  if (fp==NULL) {
    cerr << "Failed to open " << fname << endl;
    return NULL;
  }
  
  //read header
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec    = bcf_init();

  if (knownIndex > bcf_hdr_nsamples(hdr)) {
    cerr << endl << "You asked for sample number: " << knownIndex << ", but your VCF file doesn't have that many samples!" << endl << endl;
    return NULL;
  }

  
  string chrom;
  int curRid=-1;
  int a1, a2;
  int32_t ngt, *gt_arr=NULL,ngt_arr=0;
  
  bool hasKnown=knownIndex>0;
  if (hasKnown) {
    knownIndex=(knownIndex-1)*2; // genotypes stored as flattened 2D arrays (TODO, double check w/ multisample BCF)
  }
  
  //save for each vcf record
  while ( bcf_read(fp, hdr, rec)>=0 ) {
      //unpack for read REF,ALT,INFO,etc 
    bcf_unpack(rec, BCF_UN_STR);
    bcf_unpack(rec, BCF_UN_INFO);

    if (bcf_is_snp(rec) && rec->n_allele == 2) {
      Locus loc;

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
      loc.allele1 = rec->d.allele[0][0]; // by definition, must be 1 char long
      loc.allele2 = rec->d.allele[1][0];


      // extract out the genotype
      if (hasKnown) {
	ngt=bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);//FINISH 
	if (ngt <= knownIndex+1)
	  continue;
	
	a1 = bcf_gt_allele(gt_arr[knownIndex]);
	a2 = bcf_gt_allele(gt_arr[knownIndex+1]);

	// skip no-calls in the "known"
	if (a1 < 0 || a2 < 0)
	  continue;

	if (a1+a2==0) {
	  loc.genotypecall=AA;
	} else if (a1 + a2 == 1) {
	  loc.genotypecall=AB;
	} else if (a1 + a2 == 2) {
	  loc.genotypecall=BB;
	} else {
	  cerr << "Parse error for snp at position: " << loc.region << endl;
	  continue;
	}
	
      } else {
	loc.genotypecall=NOCALL;
      }
      
      loci.push_back(loc);
    }

  }
    
  bcf_destroy(rec);
  // bcf_hdr_destroy(hdr); // must be freed later. TODO!
  if (gt_arr!=NULL)
    free(gt_arr);
  
  int ret;
  if ( (ret=hts_close(fp)) ) {
    cerr << "hts_close(" << fname << "): non-zero status" << endl;
    return NULL;
  }
  return hdr;
}

bool
readBed(char *filename, vector<Locus> &loci) {
  ifstream bedFile;
  string line;
  
  bedFile.open(filename);
  if (bedFile.fail()) 
    return false;

  bool ret=true;
  bool printWarn=false;
  
  while (getline(bedFile, line)) {
    stringstream ss(line);
    vector<std::string> records;
    string tmp;
    while(getline(ss, tmp, '\t')) {
      records.push_back(tmp);
    }
    if (records.size() >0 && records.size() != 5) {
      cerr << "Error on line: " << endl << line << endl;
      ret=false;
      break;
    }
		
    int startPosition = atoi(records[1].c_str()) + 1;
    int stopPosition=atoi(records[2].c_str());
    if (startPosition==1 || startPosition>stopPosition) {
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
    loc.allele1=records[3][0];
    loc.allele2=records[4][0];
    loc.genotypecall=NOCALL;
    
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
writeCounts(const char *outfile, vector<Locus> *loci, BaseCounter *counts, double mf, const Options *opt) {

  ofstream fh;
  fh.open(outfile);
  int i;
  double error= pow(10, opt->minBaseQuality/-10.0);

  int AAindexes[N_GENOS];
  int ABindexes[N_GENOS];
  int BBindexes[N_GENOS];
  int *idx;
  int nIndexes=0;
  
  if (opt->knowns) {
    char type=MINOR;
    if (mf >= 0.5)  // major contributor is the known (double check the polarity!!)
      type=MAJOR;
    
    nIndexes = getGenoIterators("AA", type, AAindexes);
    getGenoIterators("AB", type, ABindexes);
    getGenoIterators("BB", type, BBindexes);
    
  }
  
  // convert genotypes + mixture fractions to allele weights
  double aweights[N_GENOS];
  genotypesToAlleleWeights( mf, aweights);

  double genolikes[N_GENOS];
  
  if (fh.is_open()) {
    fh << "Region\tA\tB\tMixFrac\tRefCount\tAltCount\tOtherCount\tBadCount";

    if (opt->knowns) {
      for (i =0; i < N_GENOS; ++i) {
	fh << "\t" << POSSIBLE_GENOS[i];
      }
    } else {
      fh << "\tKnown\tAA\tAB\tBB";
    }
    
    fh << endl;
    
    for (std::vector<Locus>::iterator it = loci->begin(); it != loci->end(); ++it, ++counts) {

      if (counts->badCount ==SKIP_SNP)
	continue;
      else if (opt->knowns && it->genotypecall == NOCALL)
	continue;
      
      fh << it->region << "\t" << it->allele1 << "\t" << it->allele2 << "\t" << mf <<
	"\t"  << counts->refCount <<
	"\t" << counts->altCount <<
	"\t" << counts->otherCount <<
	"\t" << counts->badCount;

      computeLikesWithM(counts->refCount, counts->altCount, aweights, error, genolikes);
      if (opt->knowns) {
	for (i =0; i < N_GENOS; ++i) {
	  fh << "\t" << genolikes[i];
	}
	
      } else {

	// get the genotype call that is known, and print it
	fh << "\t" << KNOWN_GENOS[ (int) it->genotypecall ];

	// and given what we know, grab the appropriate set of 3 genotypes of the joint likelihoods, 
	idx = AAindexes;
	if (it->genotypecall == AB)
	  idx = ABindexes;
	else if (it->genotypecall == BB)
	  idx = BBindexes;

	// and print them!
	for (i=0; i < nIndexes; ++i, ++idx) {
	  fh << "\t" << genolikes[ *idx ];
	}
	
      }
      
      fh << endl;
    }
    fh.close();
  } else {
    cerr << "Failed to open the counts file" << outfile << endl;
  }

}
