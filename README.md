# genTblastn

## Overview
This is a strategy to convert blast hits into gene models. This takes advantage of the speed of blasting. Blast hits are then parsed to give a single command to genewise, and the gff output is collected into a single file.

## Instructions
1) Blast your favorite protein set against your genome, and make use of the multithreading power of blast.

  `tblastn -query proteins.fa -db target_genome.fa -num_threads 16 -outfmt 6 > prots_vs_genome.tab`

2) Run blast2gff.py with the same three files, and output the commands automatically.

  `blast2gff.py -q proteins.fa -d target_genome.fa -b prots_vs_genome.tab -C`

This will automatically generate the command to run in parallel, which will be printed to stderr for convenience.

3) Run that command to call genewise in parallel. It will look something like this.

  `parallel --gnu -a genewise_commands_171710.sh -j 16 --joblog target_genome_genewise.log --halt 1`

4) If parallel is not installed, or cannot be used, then run command 2 without the `-C` option, to run each genewise command sequentially. For a normal sized genome with 20k proteins, this would take about 6 hours.

## History
For *ab initio* gene prediction of a new genome, it is often useful to confirm (or even just find) genes that may not be expressed (thus have no mRNA evidence) but are highly similar to known genes in other genomes. These might include developmentally restricted genes (hopefully most of them), paralogs that do not map correctly, or pseudogenes.

The goal is to convert tblastn results into gene models for use in evidenceModeler or similar evidence collection software for generation of gene models. Blast is very fast, which is why it is useful. I had looked into several other programs which generate .gff format gene models from proteins mapping onto genomes. This included [exonerate](https://www.ebi.ac.uk/~guy/exonerate/) and [genewise](http://dendrome.ucdavis.edu/resources/tooldocs/wise2/doc_wise2.html). Both of these had problems. Exonerate is unbelievably slow, and genewisedb maxed out my memory (32GB) very quickly when searching for a set of genes across the whole genome (15k prots vs. 150Mb genome contigs). Genewise also has a [difficult installation](http://ninebysix.blogspot.de/2012/11/quick-note-genewise-and-glib.html), in that it probably will not compile out of the box on linux. On ubuntu, it can be installed without downloading the source with `sudo apt-get install wise wise-doc`.

My original idea was to fork from the open source blast-to-gff program- [genBlastg](http://genome.sfu.ca/genblast/download.html), which was based on blastall or wublast. Neither of those are typically used (or even updated), so it should instead be compatible with [NCBI blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (currently v2.2.30).

As far as I could tell from looking at the code, the program creates a blast output in normal text format (the default) and then parses the text to create a "report" that vaguely resembles the tabular blast output. That is then parsed to make the genes. I had considered configuring the program to work with blast+ (which was not so hard) but the parsing of the text was so messy that I had no interest in continuing by editing existing code. The program should have taken blast output format 6 as a raw input (or generated it on the fly), but this will probably never be implemented (at least not by me)

Both versions of genBlastG that I tried (1.38, and 1.39 compiled from source) did not work; they both hit an error and died, maybe halfway through.

## Misc
Original reference for genBlastG can be found [here](http://bioinformatics.oxfordjournals.org/content/27/15/2141.full)
