# Changelog

All notable changes to Demixtify are recorded in this document.

## Post-publication
Add the ability to filter (drop reads) based on their (inferred) template length
Since I'm using up the alphabet a bit too quickly, use: <br>
	
```-z minTlenSize```
```-Z maxTlenSize```

Implementation note:
We assume that single-end reads are the result of a read merger.
Notable, single-end reads have a "0" reported tlen.
As such, we take the read length as a lower-bound on the template length
(and enforce this).

	
Note this feature helps if your samples in the mixture have different degradation profiles.
