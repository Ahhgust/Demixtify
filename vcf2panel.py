#!/usr/bin/env python3
#
# This computes the Weir and Hill (2002) estimator of Fst
# The formulas come from: Estimating and interpreting FST: The impact of rare variants
# (the SOM)
# this only works for biallelic sites!

# page 12 of the SOM

import sys
import os

def get_pbar(n1, n2, p1, p2):
    return ((n1*p1) + (n2*p2)) / (n1 + n2)


# only print snps for which the max fst (across all pop-pairs is <=:
MAX_FST=99#0.05 # disabled for his purposes...

minss= 20

PRECISION=6 # used in the round

for line in sys.stdin:
    if len(line) == 0 or line[0] == '#':
        if line.startswith("#CHROM"):
            print('##INFO=<ID=maxFst,Number=A,Type=Float,Description="Maximum pairwise FST across all population pairs">')
            print('##INFO=<ID=meanFst,Number=A,Type=Float,Description="Mean pairwise FST across all population pairs (average of ratios)">')
            print('##INFO=<ID=overallFst,Number=A,Type=Float,Description="Overall FST across all population pairs (ratio of averages)">')
            sp = line.rstrip().split("\t")
            print( "\t".join( sp[0:8]))
        else:
            print(line, end="")
        continue


    s = line.rstrip().split("\t")

    if len(s[3]) > 1 or len(s[4]) > 1:
        continue

    annos = s[7].split(";")
    
    counts = {}
    freqs = {}

    globalfreq = None
    for a in annos:
        
        if a.startswith("AN_"):
            tag, count = a.split("=")
            if len(tag)==6:
                count = int(count)
                if count >= minss:
                    counts[ tag[ 3: ] ] = count # AF_nfe -> nfe

        elif a.startswith("AF"):
            tag, count = a.split("=")
            if len(tag)==6:
                freqs[ tag[ 3: ] ] = float(count) # AF_nfe -> nfe
            elif tag == "AF":
                globalfreq = float(count)

    pops = [p for p in counts.keys() if p != "oth" and p != "raw"]
    npops = len(pops)
    
    if globalfreq is None:
        continue

        
    maxFst=-1
    npairs=0
    sumFst=0.

    numSum=0.
    denomSum=0.

    
    for i in range(npops):
        pop1 = pops[i]

        if pop1 not in freqs:
            continue

        count1 = counts[pop1]
        freq1 = freqs[pop1]
        
        for j in range(i+1, npops):
            pop2 = pops[j]

            if pop2 not in freqs:
                continue

            npairs+=1
            
            count2 = counts[pop2]
            freq2 = freqs[pop2]

            pbar = get_pbar(count1, count2, freq1, freq2)

            M = count1/count2

            numer = (1+M*M)*( (freq1-freq2)*(freq1-freq2) )
            denom = ((1+M)*(1+M)) * (pbar*(1-pbar))

            numSum += numer    
            denomSum += denom

            Fst = -1
            if denom > 0:
                Fst = numer/denom
                if Fst > 1:
                    Fst=1

                

            
            if npairs==1 or Fst > maxFst:
                maxFst=Fst
                
            sumFst += Fst


    # report ratio of averages (really sums, but that it is equivalent)
    # the max
    # and the average of ratios
    if maxFst < MAX_FST and npairs>9:      
        print("\t".join(s[0:8]), ";overallFst=" , round(numSum/denomSum, PRECISION)  , ";maxFst=", round(maxFst, PRECISION), ";meanFst=", round(sumFst/npairs, PRECISION), sep="")
    
    
