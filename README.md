# genTblastn
fork of the open source blast-to-gff program- [genBlastg](http://genome.sfu.ca/genblast/download.html)

## Overview
This is a rework of the original code, which was based on blastall or wublast. Neither of those are typically used (or even updated), so it should instead be compatible with [NCBI blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (currently v2.2.30)

Current goals are:

1. Configure to work with blast+

2. Allow skipping of the blast step with pre-blasted -outfmt 6 files 

## Misc
Original reference for genBlastG can be found [here](http://bioinformatics.oxfordjournals.org/content/27/15/2141.full)
