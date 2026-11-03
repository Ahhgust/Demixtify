# Changelog

All notable changes to Demixtify are recorded in this document.

## Post-publication
### 11/03
Small quality of life improvements
-- fixed htslib warnings (for AD and AF fields). They were "correct", but atypical.
-- AF and FST tags are no longer needed when you're DEMIX-ing (only applies when DEMIX is run in isolation)

	
### 10/31
Add the ability to filter (drop reads) based on their (inferred) template length <br>
Since I'm using up the alphabet (ie, all the flags) a bit too quickly, use: <br>
	
```-z minTlenSize```
```-Z maxTlenSize```

Implementation note:
We assume that single-end reads are the result of a read merger.
Notable, single-end reads have a "0" reported tlen.
As such, we take the read length as a lower-bound on the template length
(and enforce this).

	
Note this feature helps if your samples in the mixture have different degradation profiles.
