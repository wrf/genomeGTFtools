# genomeGTFtools

## Overview
These are some scripts to convert various features and annotations into a GFF-like file for use in genome browsers. There are also general purpose tools for genome annotation.

GFF (generic feature format) is a [tab-delimited pseudoformat](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) for annotating genomes, consisting of 8 well-defined columns, and a free-text 9th column (with some "guidelines" of what to include).

| Column | Name | Description |
|--------|------|-------------|
| 1      | seqid | The sequence ID, usually this is the contig, scaffold or chromosome (e.g. `ctg0123`, `scaffold_99`, `chrX`), and NOT the name of the gene or feature (which goes in column 9) |
| 2      | source | Usually the program that generated the data (e.g. `blastp`, `AUGUSTUS`) |
| 3      | type | Type of feature (e.g. `gene`, `exon`), typically using [standarized terms](https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/so.obo) |
| 4      | start | First base in the feature |
| 5      | end   | Last base the feature |
| 6      | score | Float (of an arbitrary scale), could be something like coverage (say 0-1000), percent identity (0.0-100.0) |
| 7      | strand | Which DNA strand, as forward (+), reverse (-), unstranded (.) or unknown (?) |
| 8      | phase | Phase of the coding sequence as 0, 1, or 2 (i.e. whether the exon ends mid-codon), only applies to `CDS` features |
| 9      | attributes | All other information, as a string of pairs of `key=value;`, though this is variable depending on version or program. Most problems with GFF are problems with parsing information from this column. |

Many of these tools were used in [our analysis of the genome of the sponge *Tethya wilhelma*](https://bitbucket.org/molpalmuc/sponge-oxygen). Please cite the paper: [Mills, DB. et al (2018) The last common ancestor of animals lacked the HIF pathway and respired in low-oxygen environments. *eLife* 7:e31176.](https://doi.org/10.7554/eLife.31176)

### Jump to: ###
* [pfam2gff.py](https://github.com/wrf/genomeGTFtools#pfam2gff) PFAM domains of proteins/coding sequences made into a GFF
* [blast2gff.py](https://github.com/wrf/genomeGTFtools#blast2gff) blast hits to GFF, generally indicating exons or conserved domains
* [microsynteny.py](https://github.com/wrf/genomeGTFtools#microsynteny) using two GFF files of two different genomes and blast hits, determine blocks of conserved gene order

## number_contigs_by_length
By convention, the longest chromosomes are numbered first. This naturally applies to scaffolds as well. Contigs/scaffolds can be renumbered and reordered with `number_contigs_by_length.py` script. Use the option `-c` to specify an additional output file of the conversion vector, that can be used to rename the scaffold column in any GFF file with the `rename_gtf_contigs.py` script.

## pfam2gff
This has two modes: one will convert the "tabular" hmmscan output (generated using PFAM-A (`Pfam-A.hmm`), which can be found in [the FTP section of PFAM](http://pfam.xfam.org/) as the hmm database, or direct download from `ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`) into a protein GFF with domains at the protein positions.

  `hmmscan --domE 0.1 --cpu 4 --domtblout stringtie.pfam.tab ~/PfamScan/data/Pfam-A.hmm stringtie_transdecoder_prots.fasta > stringtie.pfam.log`

  `pfam2gff.py -i stringtie.pfam.tab > stringtie.pfam.gff`

### For genomic coordinates ###
The other output will convert the domain positions into genomic coordinates for use in genome browsers, so individual domains can be viewed spanning exons. Run `hmmscan` as above, then use the `-g` option to include genomic coordinates. Use `-T` for presets for [TransDecoder genome GFF](https://github.com/TransDecoder/TransDecoder/wiki) file.

  `pfam2gff.py -g stringtie_transdecoder.gff -i stringtie.pfam.tab -T > stringtie_transdecoder_pfam_domains.gff`

![renilla_pfam_example.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/renilla_pfam_example.png)

For `AUGUSTUS` proteins (using [extract_features.py](https://bitbucket.org/wrf/sequences/src/master/extract_features.py) or translated nucleotides), this would be run as:

  `hmmscan --domE 0.1 --cpu 4 --domtblout renilla_test_prots.pfam.tab ~/db/Pfam-A.hmm renilla_test_prots.fasta > /dev/null`

Then run with the `AUGUSTUS` GFF (ensure this is GFF format and not GTF), instructing to use CDS features as exons with `-x`. Depending on version, the `AUGUSTUS` output might appear as below, where ID is not explicitly given in the attributes of the gene or transcripts (it just says `g1`), but is given for the introns or CDS. This makes it difficult to parse in a standarized way.

```
jcf7180000021585	AUGUSTUS	gene	921	9763	0.03	+	.	g1
jcf7180000021585	AUGUSTUS	transcript	921	9763	0.03	+	.	g1.t1
jcf7180000021585	AUGUSTUS	intron	1118	2513	0.83	+	.	transcript_id "g1.t1"; gene_id "g1";
jcf7180000021585	AUGUSTUS	CDS	921	1117	0.57	+	0	transcript_id "g1.t1"; gene_id "g1";
jcf7180000021585	AUGUSTUS	CDS	2514	2647	1	+	1	transcript_id "g1.t1"; gene_id "g1";
```

This should be changed to appear as GFF format. Introns are ignored anyway (thus can be removed) and CDS features will be treated as exons (with option `-x`).

```
jcf7180000021585	AUGUSTUS	gene	921	9763	0.03	+	.	ID=g1;Name=g1
jcf7180000021585	AUGUSTUS	mRNA	921	9763	0.03	+	.	ID=g1.t1;Parent=g1;Name=g1.t1
jcf7180000021585	AUGUSTUS	CDS	921	1117	0.57	+	0	Parent=g1.t1
jcf7180000021585	AUGUSTUS	CDS	2514	2647	1	+	1	Parent=g1.t1
```

This can be run with `pfam2gff.py`, specifying the GFF file with `-g` (this will also make it print a GFF of genomic coordinates).

  `pfam2gff.py -g aug_nocomments_cds.gff -e 1e-3 -i aug_prots.pfam.tab -x > aug_prots_pfam.gff`

Produces the following example output. The results are not filtered by position, so it is evident that 3 different von Willebrand domains are found, though group 1 is clearly the strongest match, as indicated by the score column.

```
jcf7180000021585	hmmscan	PFAM	1011	1117	113.3	+	.	ID=g1.t1.VWA.1;Name=PF00092.VWA.von_Willebrand_factor_type_A_domain
jcf7180000021585	hmmscan	PFAM	2514	2585	113.3	+	.	ID=g1.t1.VWA.1;Name=PF00092.VWA.von_Willebrand_factor_type_A_domain
jcf7180000021585	hmmscan	PFAM	1014	1117	55.5	+	.	ID=g1.t1.VWA_2.1;Name=PF13519.VWA_2.von_Willebrand_factor_type_A_domain
jcf7180000021585	hmmscan	PFAM	2514	2526	55.5	+	.	ID=g1.t1.VWA_2.1;Name=PF13519.VWA_2.von_Willebrand_factor_type_A_domain
jcf7180000021585	hmmscan	PFAM	1011	1117	35.1	+	.	ID=g1.t1.VWA_3.1;Name=PF13768.VWA_3.von_Willebrand_factor_type_A_domain
jcf7180000021585	hmmscan	PFAM	2514	2566	35.1	+	.	ID=g1.t1.VWA_3.1;Name=PF13768.VWA_3.von_Willebrand_factor_type_A_domain
```

## pfamgff2clans
Convert a PFAM protein GFF (above) to the PFAM clans, and remove some redundant hits, essentially just changing the names of the domains and merging duplicates. This is needed for the `pfampipeline.py` script. This script requires the [clan links file, called Pfam-A.clans.tsv](http://pfam.xfam.org/).

**BioPython is required for this step, to get the length of each sequence.**

## pfampipeline
**Requires BioPython, HMMER, PFAM-A and SignalP**
Script to generate graph of domains for a FASTA file of proteins. Both of the above scripts (`pfam2gff.py` and `pfamgff2clans.py`) are automatically called, followed by [SignalP](http://www.cbs.dtu.dk/services/SignalP/) and the R script `draw_protein_gtf.R`. This is called as a command on a protein file, for example on the test dataset of [nidogen](http://www.uniprot.org/uniprot/P14543) proteins:

`pfampipeline.py test_data/nidogen_full_prots.fasta`

Several output files are automatically generated, including the domain assignments in a GFF-like format, and a PDF of the domains, where a black line indicates the protein length, black box is the signal peptide, and colors correspond to different domains. **Note that some domains may overlap, due to bad calls in the HMM. Domains are by definition non-overlapping, thus these must be removed. Automatic removal may be implemented in the future.**

![nidogen_full_prots.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/nidogen_full_prots.png)

To view exon structure from a GFF file on a 3D protein structure, see instructions at my [PDBcolor repo](https://github.com/wrf/pdbcolor#gene-structure).

## blast2gff
This was a strategy to convert blast hits into gene models. The direction of the blast hit and the grouping of blast hits in the same region is most indicative of a gene (though possibly pseudogenes as well). In general, blasting all human proteins against the target genome can find many proteins even in distantly related organisms. Repeated domains or very common domains (like ATP binding for kinases) will show up all over the place, so limiting the `-max_target_seqs` is advisable.

1) Blast a protein set against the genome, and make use of the multithreading power of blast.

  `tblastn -query proteins.fa -db target_genome.fa -num_threads 16 -outfmt 6 > prots_vs_genome.tab`

2) Run blast2gff.py on the output file. This will reformat the tabular blast hits into gff3 alignment style. Simple filtering options can be applied with `-e` `-s` and `-F`. If the queries were from SwissProt, use the `-S` option to correctly format the SwissProt fasta headers for the output.

  `blast2gff.py -b prots_vs_genome.tab > prots_vs_genome.gff3`

## blast2genomegff
As above for the domains in `pfam2gff.py`, entire blast hits to transcripts or proteins can be printed as GFF using the `blast2genomegff.py` script, where the blast hits will be shown spanning multiple exons as a single feature. This might help to identify erroneously fused genes. Using transcripts from a *de novo* transcriptome [Trinity](http://trinityrnaseq.github.io/), or genome guided [StringTie](http://ccb.jhu.edu/software/stringtie/), convert blastx protein matches into genomic coordinates. For Trinity transcripts, the coordinates on the genome need to be determined by mapping the transcripts to the genome. This can be done with [GMAP](http://research-pub.gene.com/gmap/). Generically, this would be run as:

  `blastx -query transcripts.fasta -db uniprot_sprot.fasta -num_threads 16 -outfmt 6 -max_target_seqs 5 > transcripts_sprot.tab`
  
  `blast2genomegff.py -b transcripts_sprot.tab -d uniprot_sprot.fasta -g transcripts.gtf > transcripts_sprot.genome.gff`

#### Options are: ####
   * `-p` : program, by default is `blastx`, but change to `blastp` if proteins were used. The correct blast program is needed to calculate the intervals correctly, as blastp results need to be converted to nucleotide coordinates.
   * `-x` : if starting from proteins, and the `CDS` features are available (such as from `AUGUSTUS`), use CDS features from the GFF instead of exons (that is, exons are not specified at all). If exons are specified as well and are mostly identical to CDS, then also use the option `--skip-exons`.
   * `-D` : delimiter for spliting names in the tabular blast output. For example, if names were `gene1.CDS`, the option `-D .` would be used to split at `gene1`.
   * `-F` : as above for `-D`, but to split the IDs in the GFF. If the query column in blast (first column) and the ID in the GFF do not match, and there is no output, then this or `-D` may fix the problem.
   * `-c` : coverage cutoff, remove queries where the hit is under 0.1 of the subject length
   * `-e` : E-value cutoff, by default is 1e-3
   * `-s` : bitscore/length cutoff, remove hits with bitscore/length of under 0.1, that is, remove very distant matches. Set higher for more closely related species (0.3) or lower for distance species (0.05).
   * `-G` : ignore `gene` features, or other high-level features like `mRNA` or `transcript`, and instead extract the ID directly from each exon. This might be more convenient to use if exons are given unique IDs in the GFF, like `gene1.t1.exon1`. This is typically used with `-F`, as `-G -F "."`

### starting from StringTie transcripts
[StringTie](https://ccb.jhu.edu/software/stringtie/) transcripts can be converted to fasta using the script `cufflinks_gtf_genome_to_cdna_fasta.pl` (packaged with [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)). These are used as input for `blastx`. Note that with `blastx`, some can hit antisense, which suggests there is a protein on the antisense strand, or possibly there is an erroneous fusion of two adjacent genes.

1) Blastx the transcriptome against a protein set, such as [Swissprot](https://www.uniprot.org/downloads), or perhaps gene models of a related organism.

  `blastx -query stringtie.fasta -db uniprot_sprot.fasta -max_target_seqs 10 -evalue 1e-3 > stringtie_vs_swissprot_blastx.tab`

2) Convert to genomic coordinates, so individual protein hits can be seen spanning exons. The genes `-g` can be defined using the GTF output of StringTie.

  `blast2genomegff.py -b stringtie_vs_swissprot_blastx.tab -g stringtie.gtf -d uniprot_sprot.fasta -S > stringtie_vs_swissprot_blastx.gff`

### starting from AUGUSTUS proteins or CDS
When running `AUGUSTUS`, include the options `--protein=on --cds=on --gff=on`. As above, the proteins themselves (and the CDS nucleotides) can be extracted with the [extract_features.py script](https://bitbucket.org/wrf/sequences/src/master/extract_features.py). Then run either `blastp` or `blastx`, for proteins or nucleotide CDS, respectively.

  `blastp -query renilla_test_prots.fasta -db ~/db/human_uniprot.fasta -outfmt 6 -evalue 1e-4 -max_target_seqs 5 > renilla_vs_human_blastp_1e-4.tab`

By default, `AUGUSTUS` does not report exon features (only `intron` and `CDS`). 

  `blast2genomegff.py -b renilla_vs_human_blastp_1e-4.tab -g renilla_test_prots.gff -d ~/db/human_uniprot.fasta -S -x > renilla_vs_human_blastp_1e-4.gff`

If proteins were used, set `-p` to `blastp`. If nucleotide CDS was used, the default `-p` is `blastx`.

![renilla_blast_v_human_example.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/renilla_blast_v_human_example.png)

## microsynteny
Blocks of colinear genes between two species can be identified using `blast` and GFF files of the gene positions for each species.

#### Required terms are: ####
   * `-b` tabular blast or diamond results, of the query proteins or genes against some database
   * `-q` GFF file of the query genes. **The results are sensitive to presence of multiple features at the same locus. It is advisable to include only features that would match the blast results, either gene or mRNA.** i.e. `grep gene annotation.gff > genes_features_only.gff`
   * `-d` GFF of the genes of the target genome

#### Common options are: ####
   * `-m` minimum length of a block, 3 is usually sufficient for true synteny, as determined by randomized gene order (with `-R`)
   * `-z` max allowed distance between genes. This should be roughly the upper limit of intergenic distances in the genome, meaning 99% of genes are closer than `-z` (see [here for an example](https://github.com/wrf/misc-analyses/tree/master/intron_evolution)).
   * `-R` randomize gene order of the query, to estimate false discovery rate
   * `--make-gff` produce a GFF output, instead of gene-by-gene information of the synteny blocks

1) Blastx or blastp of the transcriptome (or translated CDS) against a database of proteins from the target species. Use the tabular output `-outfmt 6`. A maximum number of sequences does not need to be set, since the objective is to find homology, and this is not assumed from blast similarity. For example, here I am using the genomes of the corals [Acropora digitifera](http://marinegenomics.oist.jp/coral/viewer/info?project_id=3) and [Styllophora pistillata](http://spis.reefgenomics.org/).

  `blastp -query adi_aug101220_pasa_prot.fa -db Spis.genome.annotation.pep.longest.fa -outfmt 6 -evalue 1e-2 > acropora_vs_styllophora_blastp.tab`

2) Determine the blocks with the `microsynteny.py` script. Various delimiter parameters may need to be set depending on the version of GFF for the annotations and the names of the proteins. In this case, as the subject proteins contain the same names as the `ID` field in the GFF, the delimiter `-D` must be set to any alternate symbol, here using `|`.

  `microsynteny.py -b acropora_vs_styllophora_blastp.tab -q test_data/adi_aug101220_pasa_mrna_t1_only.gff -d test_data/Spis.genome.annotation.mrna_only.gff -D "|" > test_data/acropora_vs_styllophora_microsynteny.tab`

The standard output is a tab-delimited text file of 11 columns. The colinear block is evident from the systematic numbering of the two genomes, as the five *A. digitifera* genes (numbered 23536-23540) correspond to the five *S. pistillata* genes (Spis14158-14162). The strand direction is the same for each gene as well.

```
query-scaf	sub-scaf	block-id	query-gene	q-start	q-end	q-strand	sub-gene	s-start	s-end	s-strand
scaf16952	Spis.scaffold280|size416857	blk-2	aug_v2a.23536.t1	23026	25506	+	Spis14158	255619	257861	+
scaf16952	Spis.scaffold280|size416857	blk-2	aug_v2a.23537.t1	28396	43517	+	Spis14159	264153	287169	+
scaf16952	Spis.scaffold280|size416857	blk-2	aug_v2a.23538.t1	53236	67680	+	Spis14160	294209	310911	+
scaf16952	Spis.scaffold280|size416857	blk-2	aug_v2a.23539.t1	68035	73535	-	Spis14161	321061	322080	-
scaf16952	Spis.scaffold280|size416857	blk-2	aug_v2a.23540.t1	77256	88420	+	Spis14162	340030	347635	+
```

There are two main caveats to the data generated by this program. First, the reported block length represents the most genes on *either query or target* scaffolds, meaning if genes are erronously fused or split on one of the two species, this could, for instance, generate a block of 1 gene in the query and 3 genes in the target. Below is an example of a block where gene `aug_v2a.19447.t1` blasts to five genes in a row in *S. pistillata*. From these data alone, it cannot be determined if the *A. digitifera* gene is a false fusion, or the *S. pistillata* genes are falsely split.

```
scaf11463	Spis.scaffold21|size1341190	blk-19	aug_v2a.19446.t1	188491	189125	-	Spis2568	665897	666978	-
scaf11463	Spis.scaffold21|size1341190	blk-19	aug_v2a.19447.t1	192487	193775	-	Spis2569	668955	676058	-
scaf11463	Spis.scaffold21|size1341190	blk-19	aug_v2a.19447.t1	192487	193775	-	Spis2573.t2	740727	760040	-
scaf11463	Spis.scaffold21|size1341190	blk-19	aug_v2a.19447.t1	192487	193775	-	Spis2572	715149	725680	-
scaf11463	Spis.scaffold21|size1341190	blk-19	aug_v2a.19447.t1	192487	193775	-	Spis2571	705673	706617	-
scaf11463	Spis.scaffold21|size1341190	blk-19	aug_v2a.19447.t1	192487	193775	-	Spis2575	784361	787164	-
```

The other potential problem occurs in the case of tandem duplications. If a tandem duplication occurred in one of the two genomes, then nearly identical blocks will be identified for each duplicate, as up to 5 intervening genes are allowed (by the option `-s`). A hypothetical example is shown in the figure below. In such a case, hypothetical genes 1-3 would be matched twice (once normally, once with the intervening gene 1b), meaning the block will be counted twice.

![microsynteny_example_v1.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/microsynteny_example_v1.png)

3) Optionally, determine the frequency of random matches using the `-R` option. This will randomly reorder the positions of the query proteins, and then check for synteny. Blocks of 3-in-a-row are extremely rare, though blocks of two can be found easily by setting `-m 2`. This suggests that 3 colinear genes is usually sufficient to infer inherited gene order between the two species. Usually, long blocks are actually erroneous gene fusions in one of the genomes (as above).

  `microsynteny.py -b acropora_vs_styllophora_blastp.tab -q test_data/adi_aug101220_pasa_mrna_t1_only.gff -d test_data/Spis.genome.annotation.mrna_only.gff -D "|" -R -m 2 > test_data/acropora_vs_styllophora_microsynteny.randomized.tab`

4) A summary plot of block lengths can be produced using the associated R script `synteny_block_length_plot.R`.

  `Rscript synteny_block_length_plot.R acropora_vs_styllophora_microsynteny.tab`

![acropora_vs_styllophora_microsynteny.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/acropora_vs_styllophora_microsynteny.png)

5) Generate a GFF of the blocks, to plot in a genome browser. Use the option `--make-gff`.

  `microsynteny.py -b acropora_vs_styllophora_blastp.tab -q test_data/adi_aug101220_pasa_mrna_t1_only.gff -d test_data/Spis.genome.annotation.mrna_only.gff -D "|" > test_data/acropora_vs_styllophora_microsynteny.gff`

## scaffold_synteny
Generate a PDF of a [dot plot](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)), similar to what was done in [Srivistava 2008](https://doi.org/10.1038/nature07191) and [Simakov 2013](https://doi.org/10.1038/nature11696). This requires unidirectional blast results (not reciprocal) as duplicated blocks can be identified this way.

`scaffold_synteny.py -b hoilungia_vs_trichoplax_blastp_e-3.tab -q Hhon_BRAKER1_genes.gff3 -d Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff -f Hhon_final_contigs_unmasked.fasta -F Triad1_genomic_scaffolds.fasta --blast-query-delimiter . --blast-db-delimiter __ -l 80 -L 100 > hoilungia_vs_trichoplax_scaffold2d_points.tab`

Then, the R script is run to generate the dot plot. Synteny is clearly evident for major sections of chromosomes, indicated by the diagonals formed by the points. For instance, the longest *H. hongkongensis* contig maps completely to scaffold 7 in *Trichoplax*. This is mostly the case for the first 5 contigs in *H. hongkongensis*.

`Rscript synteny_2d_plot.R hoilungia_vs_trichoplax_scaffold2d_points.tab Hoilungia-hongkongensis Trichoplax-adhaerens`

The same can be generated for more distant species. Here two choanoflagellates are used, *Monosiga brevicollis* ([using the AUGUSTUS reannotation](https://bitbucket.org/wrf/genome-reannotations/downloads/Monbr1_augustus_v1.prots.fasta.gz)) and *Salpingoeca rosetta* ([at Ensembl](http://jul2018-protists.ensembl.org/Salpingoeca_rosetta/Info/Index)).

`blastp -query Monbr1_augustus_v1.prot_no_rename.fasta -db Salpingoeca_rosetta.Proterospongia_sp_ATCC50818.pep.all.fa -outfmt 6 -num_threads 6 -evalue 1e-3 -max_target_seqs 100 > monbr1_vs_srosetta_blastp.tab`

`scaffold_synteny.py -b monbr1_vs_srosetta_blastp.tab -q Monbr1_augustus_v1_no_comment.gff -d ~/genomes/salpingoeca_rosetta/srosetta_mrna_only_ID_renamed.gff -f Monbr1_scaffolds.fasta -F ~/genomes/salpingoeca_rosetta/Salpingoeca_rosetta.Proterospongia_sp_ATCC50818.dna.toplevel.fa.gz -l 40 -L 50 > monbr1_vs_srosetta_scaffold2d_points.tab`

`Rscript ~/git/genomeGTFtools/synteny_2d_plot.R monbr1_vs_srosetta_scaffold2d_points.tab Monosiga-brevicollis Salpingoeca-rosetta`

## extract_coordinates
This is a script to extract feature positions from a GFF to generate a figure that roughly resembles what is seen in a genome browser, to make something that has the same effect as a screenshot but looks better. Because GFFs contain indices by scaffold position, rather than gene, the script effectively converts the tabular GFF to a format where the features are the genes and exons themselves. This allows for easy downstream processing and visualization with the associated script `draw_annotation_blocks.R`.

#### Required terms are: ####
   * `-g` one or more GFF files, each roughly representing a track in a genome browser
   * `-s` the chromosome or scaffold ID
   * `-b` the beginning position to extract from that scaffold (inclusive)
   * `-e` the end position to extract from that scaffold (also inclusive)

The output of the program is a 5-column tab-delimited file, containing gene name, feature type, start, stop, and strand. Even if multiple GFFs are used, all of the transcripts or genes will end up in the same file.

For example, studying the arrangement of the Lux operon, in the bacterium [*Aliivibrio fisheri*](https://en.wikipedia.org/wiki/Aliivibrio_fischeri), the genome and GFF are downloaded [from NCBI](https://www.ncbi.nlm.nih.gov/genome/724). 

  `extract_coordinates.py -g GCF_000011805.1_ASM1180v1_genomic.gff -s NC_006841.2 -b 1041700 -e 1051700 -p > lux_locus_short_annot.tab`

## draw_annotation_blocks
A rather unsophisticated script to plot the extracted coordinates of genes or exons from `extract_coordinates.py`. By default, it only requires the output of `extract_coordinates.py` and will generate a `.pdf` of the same name. In this example, the input file `lux_locus_short_annot.tab` yields the output `lux_locus_short_annot.pdf`.

  `Rscript draw_annotation_blocks.R lux_locus_short_annot.tab`

As there is not a straightforward way to display features, or generate a figure of whatever someone would want to convey, the output is rather quirky, ultimately meaning that it is much easier to generate the `.pdf` and import it into another program like Inkscape.

* The X-axis and scaling of the features is set by the boundaries from `-b` and `-e` in `extract_coordinates.py`, meaning to make two graphs to compare loci, the span of `-b` to `-e` should be the same
* Each gene or transcript is printed on a different Y-position, of up to 30 features (this can be changed). This means that splice variants are each on their own line, as are genes in operons.
* All features are colored the same.
* The names are taken from the `ID` tag in the GFF, which may not be the same as `Name`.

The raw output is shown below. Some feature names are cut off due to the grouping of objects. These can be ungrouped and moved.

![lux_locus_short_annot.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/lux_locus_short_annot.png)

Naturally this requires some work to make it presentable, i.e. flipping the axis, coloring each gene, renaming them from the GFF ID, which did not have meaningful IDs for most of the proteins. Thus, a final version may look more like this, maintaining the color scheme used in [Dunlap 2009](http://linkinghub.elsevier.com/retrieve/pii/B9780123739445000663).

![lux_locus_nice_looking.png](https://github.com/wrf/genomeGTFtools/blob/master/test_data/lux_locus_nice_looking.png)

## repeat2gtf
From scaffolds or masked contigs, generate a feature for each long repeat of N's or n's (or any other arbitrary letter or pattern). The most obvious application is to make a track for gaps, which is the default behavior. The search is a regular expression, so could be any other simple repeat as well - CACA, CAG (glutamine repeats).

  `repeat2gtf.py scaffolds.fasta > scaffolds_gaps.gtf`

In jbrowse, I typically change a few options in the style to make this more visually useful. The color is set to black, the height is reduced to narrow bars, and the label is the score, which is the length of the gap.

```
         "style" : {
            "className" : "feature",
            "color" : "#000000",
            "height" : 5,
            "label" : "score"
         },
```

## pal2gtf
Convert palindromic repeats from the [EMBOSS program palindrome](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/palindrome.html) into GTF features. This was meant for mitochondrial genomes, but could potentially be whole nuclear genomes.

## DEPRICATED: blast2genewise
**To get gene models from blast hits, the best strategy may be to use** `blast2gff.py` **with the option** `-A` **to convert the blast hits to** [AUGUSTUS hints](http://augustus.gobics.de/binaries/README.TXT) (which are in a GFF-like format). This is then specified in the [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) run as: `--hintsfile=geneset_vs_scaffolds.gff`

  `blast2gff.py -t CDSpart -b geneset_vs_scaffolds.tab -A > geneset_vs_scaffolds.gff`

The original idea was intended to take advantage of the speed of blasting. Blast hits are then parsed to give a single command to [Genewise](http://dendrome.ucdavis.edu/resources/tooldocs/wise2/doc_wise2.html), and the gff output is collected into a single file and reformatted for modern genome browsers. This is very similar to the strategy used by the [BUSCO pipeline](https://gitlab.com/ezlab/busco), which takes blast hits and runs AUGUSTUS on the scaffold that was hit. As AUGUSTUS can use HMM-like profiles to find specific proteins that are conserved, perhaps with more complex domain structures, this might be developed further.

