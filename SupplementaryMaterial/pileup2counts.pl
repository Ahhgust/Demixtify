#!/usr/bin/env perl
#
# Written by August Woerner
# 4/04/2016
# Once upon a time I wrote this script for a job interview.
# it parses a pileup (samtools)
# and it to just count nucleotides (taking care to interpret indels correctly)
# Modified by AW; 
# This now just reports the number of As, Cs, Gs and Ts...
# (in that order)

use strict;
use warnings;
use Getopt::Long;

# PILEUP CONSTANTS
my $CHROM_COL=0;
my $POS_COL=1;
my $REFBASE_COL = 2;
my $PILE_COL = 4;
my $QUAL_COL=5;


my $MASK_CHARACTER = 'x'; # when masking out indels from the pileup-- the character I used

my $minSequenceQuality = 20; # only use (base) sequence qualities >= this
my $minCoverage =5; # inclusive; minimum coverage
my $maxCoverage =1e6; # exclusive; max coverage
my $mismatchPercentage = 0.01; 

my $help=0;
GetOptions('q=i'=> \$minSequenceQuality,
	   'p=f' => \$mismatchPercentage,
	   'm=i'=> \$minCoverage,
	   'h' => \$help,
	   'x=i' => \$maxCoverage);


if ($help) {
    die "Correct usage : samtools mpileup mpileup_command | $0\n" ,   
    "\t-q quality (the minimum base-quality when considering bases in the pile)\n",
    "\t-p mismatchPercentage (the minimum fraction of non-reference bases required for the base to be printed)\n " ,
    "\t-m minCov (the minimum acceptable coverage (inclusive))\n",
    "\t-x maxCov (the maximum acceptable coverage (exclusive))\n", 
    "\t-h (prints this warning message)\n\n";
}



my (@s, $i, $j, $len, $c, $q, $qlen);

my @EXP  = ('A', 'C', 'G', 'T');

# read stdin or ARGV
while (<>) {
    chomp;
    @s = split /\t/;
    my ($chrom, $pos, $refbase, $pile, $quals) = 
        ($s[$CHROM_COL], $s[$POS_COL], uc($s[$REFBASE_COL]), $s[$PILE_COL], $s[$QUAL_COL]);

    if (!defined $quals) {
        warn "Line number " , $. , " appears corrupt!\n$_\n";
        next;
    }

    $pile =~ s/\^.//g; # remove all mapping qualities / read starts from the pile
    $pile =~ s/\$//g; # remove the read ending markers
    # mask out insertion / deletions
    while ($pile =~ m/[+-](\d+)/g) {
        my $position = pos($pile) # right-most position of match
            - length($1); # converted to left-most position of the match (index of +/-) 
        my $maskLength =$1 + length($1) + 1; 
        #	warn $pile , ' ';
        substr($pile, $position-1, $maskLength, $MASK_CHARACTER x $maskLength); # overwrite the indel with a bunch of x's
        #	die $pile;
    }


# parse the pile string and the qual string
    $len = length $pile;
    $qlen = length $quals;
    my $totalBases=0;
    my %baseCounts = (); # contains the frequency of each base (of sufficient quality) in the pileup
    for ( $i=$j=0; $i < $len; $i++) {
        $c = substr($pile, $i, 1); 
        next if $c eq $MASK_CHARACTER; # note no j++
        $c = uc $c; # convert it to upper case
        if ($j == $qlen) {
            die "Should never happen! String too long: Error parsing the quality from line $.\n" , $_;
        }
        $q = ord( substr($quals, $j, 1)) - 33; # the phred-scale quality of base $c in integer format
        if ($q >= $minSequenceQuality) {
            if ($c eq '.' || $c eq ',') {
                $baseCounts{$refbase}++;
                $totalBases++;
            } elsif ($c =~ m/[ACGT]/) {
                $baseCounts{$c}++;
                $totalBases++;
            } 
        }
        
        $j++; 
    }
    if ($j != $qlen) {
        die "Should never happen! String too short: Error parsing the quality from line $.\n" , $_;
    }

    # check the coverage 
    if ($totalBases >= $minCoverage && 
        $totalBases > 0) {
	print "$s[0]\t", $s[1] -1 , "\t", $s[1];
	
	foreach my $b (@EXP) {
	    my $count=0;
	    if (exists $baseCounts{$b}) {
		$count = $baseCounts{$b};
	    }

	    print "\t", $count;
	}
	print "\n";
	
    }
    
}



