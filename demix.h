#ifndef DEMIX_H_
#define DEMIX_H_

#include <limits.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cstdlib>

typedef struct {
  std::string region;
  uint32_t pos;
  char allele1;
  char allele2;
  char genotypecall;
  float af; // population allele frequency; defaults to DEFAULT_AF
  float Fst; // estimator of FST; in practice, the "mean" pairwise (expectation of ratios) works the best (as opposed to typical estimator; ratio of expectatinos)
} Locus;


// the bits and pieces needed to 
// parse a VCF file
typedef struct {
  htsFile *fp;
  bcf_hdr_t *hdr;
  bcf1_t *rec;
  hts_idx_t *idx;
  const char *fname;
} VcfIterator;


typedef struct {
  std::string filename; // the name of the file being written
  std::string file_format; // w (vcf) or wb (bcf/vcf.gz)
  htsFile *fp;
  bcf_hdr_t *hdr;
} VcfWriter;

#define NOCALL 0
#define AA 1
#define AB 2
#define BB 3

// the maximum Phred-scaled quality (possible)
#define MAXQUAL 100

#define DEFAULT_AF -1
#define MIN_AF 0.0001

// 0.02 changes
// New in this version:
// Fst corrected genotype probability estiamtion
// and the expanded likeilhood function (moving from formula 5 to formula 3 of Crysup & Woerner)
#define VERSION "0.02"

// FST ala Hudson using the formulation of:
// http://www.genome.org/cgi/doi/10.1101/gr.154831.113.
// (see the SOM)
// #define FST_HUDSON(p,q) (  ( (p-q)*(p-q) )/( (p*(1.-q) + (q*(1.-p))) )  )


// the above (commented out) macro is right (numerically) for the N/D estimator of Fst, but it's
// there is a typo in the SOM (D is defined incorrectly in the text). For one thing, it's asymmetric
// (it should be p1(1-p2) + p2(1-p1), but they wrote
// p1(1-p2) + p1(1-p2). whoops!

// To avoid confusion
// let's just stick w/ the 1-Hw/Hb formulation, which is also defined in the SOM.

#define FST_HUDSON(p,q) (  1.0-( (p*(1.0-p) + q*(1.-q)) / (p*(1.0-q) + q*(1.-p)) )  )

const std::string DEFAULT_FST_TAG = "meanFst";

// -10.0/ln(10)
// useful when converting log likelihoods (LN)
// to phred-scaled likelihoods (-10*log10(like))
const double PHRED_SCALAR = -4.3429448190325175;

const double DEFAULT_MAX_FST=0.05;


const char* DEFAULT_AF_TAG="AF";

// note: KNOWN_GENOS[AB] == "AB"
std::string KNOWN_GENOS[] = {
			     "",
			     "AA",
			     "AB",
			     "BB"};

typedef struct {
  unsigned minMapQuality;
  unsigned minBaseQuality;
  unsigned maxBaseQuality;
  unsigned ngrid; // parameter sweep; how fine?
  int minReadLength; // bam spec. lens are signed ints
  int numThreads;
  int knowns; // the INDEX of the known contributor in the VCF file (1-based)
  uint32_t filter; // SAM_FLAG filter
  uint32_t include_filter; // SAM_FLAG filter
  double mixtureFraction; // can be passed in, or estimated by ML
  double downsampleFraction; //
  bool filterIndelAdjacent;
  bool help;
  bool trustQualities;
  const char *knownBcf; // known genotypes. at most 1 known contributor
  const char *lowFstBcf;  // low FST snps; used to estimate MF
  const char *referenceFasta; // needed for cram support.
  char *bamFilename;
  char *bedFilename;
  const char *outCounts;
  std::string outVCF;
  std::string fstTag; // contains the population allele frequency TAGS (eg, AF_nfe) used to estimate FST
  double fstFilt;
  const char* AFtag;
  bool parseVcf; // the "bed" file may either be .bed (simulations) or .vcf (real data)
} Options;

typedef struct {
  char b; // basecall 
  unsigned char q; // quality;
} BQ;

typedef struct {
  std::vector<unsigned char> refQuals; // not really reference; alleles of type 1 or 2. count the number of bases of sufficient quality
  std::vector<unsigned char> altQuals; // ditto for the other (alt) allele
  unsigned otherCount; // basecounts that pass filters and are neither type 1 nor type 2 
  unsigned badCount;  // the number of reads at the site that fail the QC filters
} BaseCounter;


// not strictly necessary, but for convenience, let's encode the
// 9 possible biallelic genotype combinations for a two-person mixture
// the first two chars for for person 1, the second two are for person 2
std::string POSSIBLE_GENOS[] = {
			   "AAAA",
			   "AAAB",
			   "AABB",
			   "ABAA",
			   "ABAB",
			   "ABBB",
			   "BBAA",
			   "BBAB",
			   "BBBB"
};


// the "likelihood"  (PL-scaled) for a known genotype (e.g., 0,255,255 for a homozygous reference call)
#define KNOWN_PL 255

#define N_GENOS 9
#define N_MARGINALS 6
#define SKIP_SNP UINT_MAX
#define IGNORE_SNP (UINT_MAX-1)

// gives indexes in POSSIBLE_GENOS consistent w/ some known genotype (see below)
int getGenoIterators(const std::string &known, char type, int *itr);


// in terms of getGenoIterators, the type considered
#define EITHER 0
#define MAJOR 1
#define MINOR 2
#define ALL   3

#define GENOPROB(c, a)       ( c[0]=='A' && c[1]=='A' ? (1-a)*(1-a)             : (c[0]=='B' && c[1]=='B' ? (a)*(a)             : 2*(1-a)*a))

// theta-corrected
#define GENOPROBFST(c, a, f) ( c[0]=='A' && c[1]=='A' ? ((1-a)*(1-a)+(1-a)*a*f) : (c[0]=='B' && c[1]=='B' ? ((a)*(a)+(1-a)*a*f) : 2*(1-a)*a*(1.-f) ))

// Meat and potatoes!
// this takes in a single locus and fills out a BaseCounter struct
// it does so by taking a pileup at the location give in loc
// while adhering to the filters in opt
bool summarizeRegion(samFile *in, bam1_t *b, sam_hdr_t *header, hts_idx_t *idx, Locus *loc, BaseCounter &count, Options &opt);



// ################# IO ROUTINES
// Used to define loci. Needs to be either a bed file
// or a v/bcf file (using htslib).
// the latter is for convenience,
// while the former (bed) is more pure (the two alleles need not specify a reference allele,
// which would be implicitly assumed by read_vcf (one ref allele and one alt allele)
// returns NULL if the file cannot be read (or other IO problems)
// return type is bcf_header_t*, but to leave the .h file simple I called it void*
// Note that *results MAY BE NULL
// if it is, then just the genotypes are extracted.
// otherwise the DP field is used to populate the BaseCounter object
void *readLowFstPanel(const char *fname, std::vector<Locus> &loci, BaseCounter *results, const char* AFtag, const char* fstTag);

bool readBed(char *filename, std::vector<Locus> &loci, std::vector<BaseCounter> &tmpResults, const Options &opt);

// writes the counts data structure to file (plain text)
// and it writes the likelihoods...
void writeCounts(const char *outfile, std::vector<Locus> *loci, BaseCounter *counts, double mf, double error, const Options *opt);

// # Helper/convenience functions
// validates a an allele (must be in the set [ACGT])
// biallelic snps only
bool validNuc(char c);

// NOTE2SELF:
// using natural logs! These must be converted to log10s before the GL is computed!
// converts 'w' per crysup and woerner (2022)
// ie, the (expected) "weight" of allele A in a biallelic (alleles A and B) two-person mixture (see POSSIBLE_GENOS)
#define IS_A(A) (A=='A')
#define BIALLELIC_WEIGHT(g,f) IS_A( (g)[0])*(f)/2.0 + IS_A((g)[1])*(f)/2.0 + IS_A( (g)[2])*(1.0-(f))/2.0 + IS_A((g)[3])*(1.0-(f))/2.0;


// formula 5 from Crysup and Woerner
// note: this applies to a site (across bases)
#define LOG_LIKE(NC, NT, w, e) (NC*log(w*(1.-e) + (1.-w)*e/3.0) + NT*log( (1.-w)*(1.-e)+w*(e/3.))  )

// part of formula 3 of Crysup and Woerner
// namely: w*(Pr(b)|A)
// note: this applies to a single BASE at a site
// takes in the weight (w) (real), the qualiy score (byte), and whether or not the base matches the proposed allele (bool)
#define LIKE_PER_NUC(w, q, m) ( m ? (w)*(1.0-QUAL_2_ERRORPROB[q]) : (w)*(QUAL_2_ERRORPROB[q]/3.0)   )

// log_sum_exp on a pair of log likelihoods (ie, gives log of sum of probabilities)
// note the 1 is for taking the exp(max-max), which happens implicitly
#define LOG_SUM_EXP_PAIR(L1, L2) ( (L1)>(L2) ? log( exp((L2)-(L1)) + 1. )+(L1) : log( exp((L1)-(L2)) + 1. )+(L2) )
#define LOG_SUM_EXP_PAIR_NAIVE(L1, L2) ( log(exp(L1)+exp(L2)))

// converts all proposed genotypes  e.g., (AAAB == AA individual 1, AB indivual 2)
// and the propsed mixture fraction (eg, mf==70% of reads come from individual 1)
// into the expected proportion of reads associated with an "A" allele
// *w must be of length 9; data written to it. w is parallel to POSSIBLE_GENOS
void genotypesToAlleleWeights(double mf, double *w);

// computes the 9 likelihoods associated with a pair of counts on A and B alleles
void computeLikesWithM(const BaseCounter *b, double *w, double e, double out[N_GENOS], double altFreq, double fst, bool trustQualities);



int csv2vec(const char* s, std::vector<std::string> *vec);
int csv2vec(const char* s, std::vector<double> *vec);

// log-sum-exp trick applied to the genotype likelihoods
// gives the log( sum (likelihoods)) when the likelihoods are in log-space
// 
inline double log_sum_exp(double *genolikes);

// same as above, but a vector of indexes (knowns, corresponding to 2-person genotypes
// which are consistent with some KNOWN genotype) can be used to compute the LL
// considering SOME of the indexes
inline double
log_sum_exp_with_knowns(double *genolikes, int *knowns, int ngenos);
// perl-inspired print to stderr and exit. extra may be NULL
void die(const char *arg0, const char *extra);

// parses the argv
// optionally fills out a Locus struct
bool parseOptions(char **argv, int argc, Options &opt, Locus &loc);

// takes in a vector of loci (loc), and a bcf file, and adds in the known genotype field.
// returns the number of genotypes recovered...
// seeks to each snp position; use this when we estimate the mixture proportion
unsigned addKnowns(std::vector<Locus> &loci, const char* knownBcf, int knownIndex);

// seeks to the chromosome; linear scan to populate records
// use this for genotyping (a lot of) snps
unsigned addKnowns(std::vector<Locus> &loci, const char* knownBcf, int knownIndex, const std::string &chr);

bool getNextVcfRec(VcfIterator &vcf, Locus &loc);
VcfIterator initVcfIterator(const char *filename, int knownIndex, bool needidx);
bool closeVcfIterator(VcfIterator &vcf);

char vcfGetGenotype(VcfIterator &vcf, Locus &loc, int knownIndex);

#endif
