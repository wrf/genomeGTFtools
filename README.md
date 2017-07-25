# genomeGTFtools

## Overview
These are some scripts to convert various features and annotations into a GFF-like file for use in genome browsers.

## pfam2gff
This has two modes: one will convert the "tabular" hmmscan output (generated using [PFAM-A](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) as the database) into a protein GFF with domains at the protein positions. 

  `hmmscan --cpu 4 --domtblout stringtie.pfam.tab ~/PfamScan/data/Pfam-A.hmm stringtie_transdecoder_prots.fasta > stringtie.pfam.log`

  `pfam2gff.py -i stringtie.pfam.tab > stringtie.pfam.gff`

The other output will convert the domain positions into genomic coordinates for use in genome browsers, so individual domains can be viewed spanning exons. Run `hmmscan` as above, then use the `-g` option to include genomic coordinates. Use `-T` for presets for TransDecoder genome GFF file.

  `pfam2gff.py -g stringtie_transdecoder.gff -i stringtie.pfam.tab -T > stringtie_transdecoder_pfam_domains.gff`

## pfamgff2clans
Convert a PFAM protein GFF (above) to the PFAM clans, and remove some redundant hits, essentially just changing the names of the domains and merging duplicates. This is needed for the pfampipeline.py script. This script requires the [clan links](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz).

## repeat2gtf
From scaffolds or masked contigs, generate a feature for each long repeat of N's or n's (or any other arbitrary letter or pattern). The most obvious application is to make a track for gaps, which is the default behavior. The search is a regular expression, so could be any other simple repeat as well - CACA, CAG (glutamine repeats).

  `repeat2gtf.py scaffolds.fasta > scaffolds_gaps.gtf`

## pal2gtf
Convert palindromic repeats from the EMBOSS program palindrome into GTF features. This was meant for mitochondrial genomes, but could potentially be whole nuclear genomes.

## blast2gff
This was a strategy to convert blast hits into gene models. The direction of the blast hit and the grouping of blast hits in the same region is most indicative of a gene (though possibly pseudogenes as well). In general, blasting all human proteins against the target genome can find many proteins even in distantly related organisms. Repeated domains or very common domains (like ATP binding for kinases) will show up all over the place, so limiting the `-max_target_seqs` is advisable.

1) Blast a protein set against the genome, and make use of the multithreading power of blast.

  `tblastn -query proteins.fa -db target_genome.fa -num_threads 16 -outfmt 6 > prots_vs_genome.tab`

2) Run blast2gff.py on the output file. This will reformat the tabular blast hits into gff3 alignment style. Simple filtering options can be applied with `-e` `-s` and `-F`. If the queries were from SwissProt, use the `-S` option to correctly format the SwissProt fasta headers for the output.

  `blast2gff.py -b prots_vs_genome.tab > prots_vs_genome.gff3`

## blast2genomegff
Instead, using transcripts from a *de novo* transcriptome [Trinity](http://trinityrnaseq.github.io/), or genome guided [StringTie](http://ccb.jhu.edu/software/stringtie/), convert blastx protein matches into genomic coordinates.

1) Blastx the transcriptome against a protein set, such as Swissprot, or perhaps gene models of a related organism.

  `blastx -query transcripts.fasta -db uniprot_sprot.fasta -num_threads 16 -outfmt 6 -max_target_seqs 5 > transcripts_sprot.tab`
  
2) Convert to genomic coordinates, so individual protein hits can be seen spanning exons.

  `blast2genomegff.py -b transcripts_sprot.tab -d uniprot_sprot.fasta -g transcripts.gtf > transcripts_sprot.genome.gff`

For Trinity transcripts, the coordinates on the genome need to be determined by mapping the transcripts to the genome. This can be done with [GMAP](http://research-pub.gene.com/gmap/).

## DEPRICATED: blast2genewise
This was intended to take advantage of the speed of blasting. Blast hits are then parsed to give a single command to Genewise, and the gff output is collected into a single file and reformatted for modern genome browsers. Given that Genewise is depricated, this will likely be changed to use AUGUSTUS in the future in order to find specific proteins that are conserved, perhaps with more complex domain structures.

For long proteins with repeated domains, the prediction will probably not work well unless the query can cover the entire gene on the genome.

Groups of blast hits in the same region define the boundary for Genewise to speed up the gene search, plus a margin on both sides. All gff outputs of Genewise are collected into a single file that is named automatically. These will be in the normal gene-mRNA-exon-CDS format for gff3 files. 

  `blast2genewise.py -q proteins.fa -d target_genome.fa -b prots_vs_genome.tab`
  
  For some genome browsers, exon features in the gff3 may clutter up the viewing window, therefore can be excluded with the `-E` flag. In most cases it is good to have them, since they can also be removed later quite easily with `grep -v exon`.

I previously had the idea to use parallel to run a bunch of Genewise commands. However, these would all have to be wrapped again since the Genewise gff format does not provide the Name or ID of each feature, so is probably incompatible with most genome browsers.

Because similar proteins or splice variants tend to produce identical gene predictions, these can be removed with the included script `removeredundantgff.py` as:

   `removeredundantgff.py -g target_genome_genewise.gff > target_genome_genewise.unique.gff`

### A lengthy explanation
For *ab initio* gene prediction of a new genome, it is often useful to confirm (or even just find) genes that may not be expressed (thus have no mRNA evidence) but are highly similar to known genes in other genomes. These might include developmentally restricted genes (hopefully most of them), paralogs that do not map correctly, or pseudogenes.

The goal was to convert tblastn results into gene models for use in evidenceModeler or similar evidence collection software for generation of gene models. Blast is very fast, which is why it is useful. I had looked into several other programs which generate .gff format gene models from proteins mapping onto genomes. This included [exonerate](https://www.ebi.ac.uk/~guy/exonerate/) and [Genewise](http://dendrome.ucdavis.edu/resources/tooldocs/wise2/doc_wise2.html). Both of these had problems. Exonerate is unbelievably slow, and genewisedb maxed out my memory (32GB) very quickly when searching for a set of genes across the whole genome (15k prots vs. 150Mb genome contigs). Genewise also has a [difficult installation](http://ninebysix.blogspot.de/2012/11/quick-note-genewise-and-glib.html), in that it probably will not compile out of the box on linux. On Ubuntu, it can be installed without downloading the source with `sudo apt-get install wise wise-doc`.

My original idea was to fork from the open source blast-to-gff program- [genBlastg](http://genome.sfu.ca/genblast/download.html), which was based on blastall or wublast. Neither of those are typically used (or even updated), so it should instead be compatible with [NCBI blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (currently v2.2.30).

As far as I could tell from looking at the code, the program creates a blast output in normal text format (the default) and then parses the text to create a "report" that vaguely resembles the tabular blast output. That is then parsed to make the genes. I had considered configuring the program to work with blast+ (which was not so hard) but the parsing of the text was so messy that I had no interest in continuing by editing existing code. The program should have taken blast output format 6 as a raw input (or generated it on the fly), but this will probably never be implemented (at least not by me)

Both versions of genBlastG that I tried (1.38, and 1.39 compiled from source) did not work; they both hit an error and died, maybe halfway through. The code appeared to be last updated in 2012, so this may not be under development anymore. The original reference for genBlastG can be found [here](http://bioinformatics.oxfordjournals.org/content/27/15/2141.full).

## Misc
This is a work in progress; there may be a citation to come once anything gets used in a real paper.
