# archived scripts #

## pal2gtf
Convert palindromic repeats from the [EMBOSS program palindrome](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/palindrome.html) into GTF features. This was meant for mitochondrial genomes, but could potentially be whole nuclear genomes.

## OBSOLETE: blast2genewise
**To get gene models from blast hits, the best strategy may be to use** `blast2gff.py` **with the option** `-A` **to convert the blast hits to** [AUGUSTUS hints](http://augustus.gobics.de/binaries/README.TXT) (which are in a GFF-like format). This is then specified in the [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) run as: `--hintsfile=geneset_vs_scaffolds.gff`

  `blast2gff.py -t CDSpart -b geneset_vs_scaffolds.tab -A > geneset_vs_scaffolds.gff`

The original idea was intended to take advantage of the speed of blasting. Blast hits are then parsed to give a single command to [Genewise](http://dendrome.ucdavis.edu/resources/tooldocs/wise2/doc_wise2.html), and the gff output is collected into a single file and reformatted for modern genome browsers. This is very similar to the strategy used by the [BUSCO pipeline](https://gitlab.com/ezlab/busco), which takes blast hits and runs AUGUSTUS on the scaffold that was hit. As AUGUSTUS can use HMM-like profiles to find specific proteins that are conserved, perhaps with more complex domain structures, this might be developed further.

Nowadays (in 2026), it is better to use long-read RNAseq plus something like `miniprot` mappings.

