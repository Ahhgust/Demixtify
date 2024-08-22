#!/usr/bin/env python3
#
# Written by August Woerner
#
# This is DNA mixture simulator.
# It strictly applies to biallelic SNPs and two-person DNA mixtures.
# The parameters of the mixture are known (specified), as are properties on the amount of
# data collected (e.g., -t, the total sample coverage).
# Read depths are assumed to follow a Poisson distribution, and
# Hardy Weinberg AND linkage equilibrium are assumed.
#
# The basic operation is as follows:
# Some allele frequency is provided (AF tag the VCF file is the default; though one or two different AF tags can be used; see -F and -E)
# and this is converted into genotype probabilities (HW)
# two genotypes are drawn (independently)
# and then number of reads observed at the site is drawn (Poisson, given the total coverage)
#         At this point we know, at a site, how many reads we have observed
# The total number of reads are then partitioned into the two individuals that constitute the mixture (binomial)
#         At this point we known how many reads belong to person 1 and person 2
# And then alleles are drawn given the true genotype (binomial)
#         At this point we known (in truth) how many A and B alleles we have sampled (without error)
# And then sequencing error is simulated (binomial, though how it is done is equivalent to a multinomial)
#         At this point we know how many A and B alleles we have measured (WITH error)
# The number of reads for alleles A and B are then summed (across the two individuals, thus mixing them
#         The above was for one person's reads; now we combine them.
#
# And we report the number of reads that support each allele (measured, As and Bs), as well as what the true genotypes are.


from numpy import random
import sys
import argparse
from signal import signal, SIGPIPE, SIG_DFL  
signal(SIGPIPE,SIG_DFL) # disable pipe errors (eg, by streaming this program to head)


def getAltFreq(info, tagname):
    sp = info.split(";")
    for tag in sp:
        dat = tag.split("=")
        if len(dat)==2:
            if dat[0]==tagname:
                return float(dat[1])

    return -1


def main(argv):

    
    parser = argparse.ArgumentParser(description="Let's make some data!\n")
    parser.add_argument('-f', '--alt_allele_frequency', dest='F', help="The (singular) alternative allele frequency modeled (assumes HWE). Default: 0.5",default=0.5, type=float)
    parser.add_argument('-F', '--vcf_alt_allele_frequency_tag', dest='TAG', help="Grabs the alternative allele frequency from the V/BCF (provide the tag; typically AF. Assumes HWE).",default="", type=str)
    parser.add_argument('-E', '--vcf_alt_allele_frequency_tag_2nd_contrib', dest='TAG2', help="Grabs the alternative allele frequency from the V/BCF for the 2nd contributor (provide the tag; defaults to -F (intrapopulation). Assumes HWE).",default="", type=str)
    parser.add_argument('-G', '--vcf_alt_allele_frequency_tag_to_filter', dest='TAGFILTER', help="Uses the alternative allele frequency from the V/BCF (to LB/UB filter), genotypes come from -F",default="", type=str)
    parser.add_argument('-L', '--vcf_alt_allele_frequency_lower_bound', dest='LB', help="Lower bound on the alt allele frequency from the V/BCF (provide the tag; typically AF. Assumes HWE). Default: 0.0",default=0.0, type=float)
    parser.add_argument('-U', '--vcf_alt_allele_frequency_upper_bound', dest='UB', help="Lower bound on the alt allele frequency from the V/BCF (provide the tag; typically AF. Assumes HWE). Default: 1.0",default=1.0, type=float)
    parser.add_argument('-w', '--wrights_fst_tag', dest='W', help="Tag for FST in the INFO field", type=str, default="maxFst")
    parser.add_argument('-q', '--quality', dest='Q', help="The (singular) Phred quality. Default: 20",default=20., type=float)
    parser.add_argument('-s', '--seed', dest='S', help="The prng seed. Default: (seeds from the clock)",default=None, type=int)
    parser.add_argument('-1', '--coverage_person1', dest='C1', help="The coverage (expected read depth) for person 1. Default 30.",default=30., type=int)
    parser.add_argument('-2', '--coverage_person2', dest='C2', help="The coverage (expected read depth) for person 2. Default 30.",default=30., type=int)
    parser.add_argument('-t', '--total_coverage', dest='T', help="The coverage (expected read depth) for the mixture (sum person1,person2). Defaults to using -1 + -2",default=-1, type=float)
    parser.add_argument('-r', '--ratio', dest='R', help="The fraction of -t/--total_coverage associated with person2",default=-1, type=float)
    parser.add_argument('-m', '--max_snps', dest='M', help="The maximum number of snps emitted. Default: All.",default=3.3e9, type=int)
    parser.add_argument('-p', '--pileup', dest='P', help="Prints a 'pileup' a la VerifyBamID. ",default=False, action='store_true')
	
    results = parser.parse_known_args(argv[1:])[0]
    args = parser.parse_known_args(argv[1:])[1]

    # this is used to mismodel the population allele frequency
    # compute genotypes wrt to TAG (eg, AF_nfe) and filter based on some other population (AF_afr)
    # if no filter population is specified, assume they are samsies
    if results.TAG != "":
        if results.TAGFILTER=="":
            results.TAGFILTER=results.TAG
        if results.TAG2=="":
            results.TAG2=results.TAG

    if len(args):
        print("Extra arguments detected...", args, sep="\n", file=sys.stderr)
        exit(1)


    #popSet = set( results.W.split(",") )

    fstTag = results.W
    
    random.seed(results.S)# optional determinism
    
    pError = 10**(results.Q/-10.0)
        
    pAlt = results.F*results.F
    pRef = (1-results.F)*(1-results.F)
    pHom = pAlt+pRef

    c1 = results.C1
    c2 = results.C2
    if results.T >=0:
        c2 = results.T*results.R
        c1 = results.T*(1.0-results.R)

    if c1 < 0 or c2 < 0:
        print("That's an error!", c1, c2, file=sys.stderr)
        return 1


    i=0
    prevpos=-1
        # takes in a VCF file from stdin
    for line in sys.stdin:
        if line[0] == '#':
            continue

        sp = line.rstrip().split("\t")
        # biallelic SNP sites only. index 3,4 refer to ref,alt allele call
        if sp[4] == '.' or len(sp[3]) != 1 or len(sp[4]) !=1:
            continue

        if i >= results.M:
            break

        fst = getAltFreq(sp[7], results.W)
        # parsing failure
        if fst < 0:
            continue
        
        if results.TAG:
            altFreq = getAltFreq(sp[7], results.TAG)
            filterFreq = getAltFreq(sp[7], results.TAGFILTER)
            if altFreq < 0: # no alternative allele frequency found...
                continue
            if filterFreq < results.LB or filterFreq > results.UB:
                continue

            
            pAlt = altFreq*altFreq
            pRef = (1-altFreq)*(1-altFreq)
            pHom = pAlt+pRef
        
        d1 = getCounts(c1, pError, pRef, pHom)
        
        if results.TAG2!=results.TAG: # both !="" and != each other -> different populations
            altFreq = getAltFreq(sp[7], results.TAG2)
            if altFreq < 0: # no alternative allele frequency found...
                continue

            pAlt = altFreq*altFreq
            pRef = (1-altFreq)*(1-altFreq)
            pHom = pAlt+pRef
        
        d2 = getCounts(c2, pError, pRef, pHom)

        pos = int(sp[1])
        if pos == prevpos:
            continue
        prevpos=pos


        nother = d1[-1] - (d1[1] + d1[2]) # total - (numAs - numBs), person1
        nother += d2[-1] - (d2[1] + d2[2]) # total - (numAs - numBs), person2
        
        if results.P: # need to print "pileup" format

            if d1[3] + d2[3] == 0:
                continue
            
            nref = d1[1] + d2[1]
            nalt = d1[2] + d2[2]

            out = []
            for j in range(nalt):
                if j % 2 == 0:
                    out.append(sp[4])
                else:
                    out.append(sp[4].lower())
			
            for j in range(nref):
                if j % 2 == 0:
                    out.append(",")
                else:
                    out.append(".")
			
            outs = "".join(out)

            if nref + nalt < c1+c2:
                otherallele=notIt(sp[3], sp[4])
                outs += otherallele * int(round( nother ))
			
            quals = chr( round(results.Q) + 33) * int(round(c1 + c2))
            print(sp[0], sp[1], sp[3], d1[3]+d2[3], outs, quals, sep="\t")
			
			


        elif results.TAG:
            tags = str(filterFreq)
            tags += "," + str(fst)
            print( sp[0], pos-1, pos, sp[3], sp[4], d1[1] + d2[1], d1[2] + d2[2], nother, tags, d1[0], d2[0], sep="\t")
            
        else:
            tags = str(results.F)
            tags += "," + str(fst)
            
            print( sp[0], pos-1, pos, sp[3], sp[4], d1[1] + d2[1], d1[2] + d2[2], nother, tags, d1[0], d2[0], sep="\t")
                
        i += 1


def notIt(a, b):
	guesses = ['A', 'C', 'G', 'T']
	
	for guess in guesses:
		if guess != a and guess != b:
			return guess

	exit(1)
	
	
def pAdd(obs1, obs2):
    return (obs1[0] + obs2[0], obs1[1] + obs2[1])

def getGeno(pRef, pHom):
    r = random.uniform()
    # 0/1/2 encoding (number of alternative alleles
    gt=1
    if r < pRef:
        gt=0
    elif r < pHom:
        gt=2
    return gt

def addError(nAlleles, perror):
    '''
    Takes some number of alleles (sampled)
    and adds in errors
    returns a tuple of length 2:
    count 1 is the number of alleles observed without error, 
    count 2 is the number of alleles observed with error
    assumes 3 possible error states with a uniform probability of observing each error state
    and assumes only 1 of these error states is associated with a "bad" (ref->alt and vice versa) error
    (other allele counts are simply dropped)
    '''
    s = random.binomial(nAlleles, perror)# how many mistakes?
    sprime=random.binomial(s, 1/3.)#of the mistaken alleles, how many changed a ref to an alt (or vice versa)
    return (nAlleles-s, sprime)

def getCounts(cov, perror, pRef, pHom):
    '''
    returns a tuple of length two, corresponding to
    an observed number of reference (index0) and alternative (index1)
    allele counts
    encorporates sequencing error and sampling error
    '''
    gt = getGeno(pRef, pHom)
    depth=random.poisson(cov)
    if depth==0:
        return (gt, 0,0, depth)

    if gt==0:
        tup = addError(depth, perror)
        return (0, tup[0], tup[1], depth)
    elif gt==2:
        tup = addError(depth, perror)
        return (2, tup[1], tup[0], depth) # flip count indexes == sampling other allele
    else:
        a1= random.binomial(depth, 0.5)# how many dad alleles did we sample?
        a2 = depth-a1 #implies number of mom alleles
        obs1=addError(a1, perror) # each allele is measured independently
        obs2=addError(a2, perror)
        both = pAdd(obs1, (obs2[1], obs2[0])) # read depths are additive
        return (1, both[0], both[1], depth)
    
if __name__ == '__main__':
    main(sys.argv)
