# genomeGTFtools

## Tips for use in [JBrowse](https://jbrowse.org/) ##
JBrowse is a genome browser that can run in a web browser. Because of the lack of consistency in GTF/GFF formats (all of which hinges on the flexibility of the final column), various genome browsers will interpret GFF files differently.

Thus, to best configure GFF files for display in JBrowse, there are several considerations.

* Parent features need to be generated for anything in the GFF spanning more than 1 line. This is for memory considerations in how the information is fetched. This means that several of the scripts will have to generate additional lines for a parent feature, rather than relying on the ID of each feature.
* Some features are given by particular programs, but relationships between them are not recognized by default by JBrowse. For instance, [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) will output a GFF with the features `gene-transcript-intron-CDS`. This needs to be converted to `gene-mRNA-CDS` in order to display correctly.
* Older GTF versions will likely need substantial format conversion to display as expected.

## An example workflow ##
Here I make use of some examples from the [genome](https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/) of the [placozoan *Hoilungia hongkongensis*](https://doi.org/10.1371/journal.pbio.2005359). I am using JBrowse version 1.8 on Ubuntu 16.04.

In preparing files for JBrowse, the default directory for the objects is `/var/www/html/jbrowse/data/`. Here, to set up subfolders to allow multiple genomes in the same JBrowse instance, the local folder is specified in the script `prepare-refseqs.pl` using the option `--out`. 

```
cd /mnt/genome_data/Hhon/
prepare-refseqs.pl --fasta Hhon_final_contigs_unmasked.fasta --out ./
```

### Adding StringTie and TransDecoder ###
Because of issues with parsing the GTF format of [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml), this needs to be converted to a GFF file with a Parent-Child feature, here simply as mRNA-exon. By default, each transcript and exon will be a standalone feature, rather than a transcript composed of exons as subfeatures.

```
stringtie_gtf_to_gff3.py Hhon_tophat2_stringtie.gtf > Hhon_tophat2_stringtie.gff3
flatfile-to-json.pl --gff Hhon_tophat2_stringtie.gff --trackType CanvasFeatures --trackLabel Stringtie --out ./
```

Several other configuration changes are used. CDS features are expected for subfeatures, but these are not specified in the GFF. JBrowse can be instructed to use exons instead, by changing the `subParts` tag.

`"subParts" : "exon",`

For [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki), the translations can be turned into a `genome-GFF` as given by the instructions. The first five lines are given below.

```
contig_001_length_2062630	transdecoder	gene	7072	9993	.	-	.	ID=stringtie.1;Name=ORF
contig_001_length_2062630	transdecoder	mRNA	7072	9993	.	-	.	ID=stringtie.1.1|m.290;Parent=stringtie.1;Name=ORF
contig_001_length_2062630	transdecoder	exon	9894	9993	.	-	.	ID=stringtie.1.1|m.290.exon1;Parent=stringtie.1.1|m.290
contig_001_length_2062630	transdecoder	exon	9565	9746	.	-	.	ID=stringtie.1.1|m.290.exon2;Parent=stringtie.1.1|m.290
contig_001_length_2062630	transdecoder	CDS	9565	9744	.	-	.	ID=cds.stringtie.1.1|m.290;Parent=stringtie.1.1|m.290
...
```

In the genome-GFF output file, the `name` tag of each protein is "ORF". This can be changed with the `change_transdecoder_names.py` script to give a more informative name. 

```
contig_001_length_2062630	transdecoder	gene	7072	9993	.	-	.	ID=stringtie.1;Name=stringtie.1
contig_001_length_2062630	transdecoder	mRNA	7072	9993	.	-	.	ID=stringtie.1.1|m.290;Parent=stringtie.1;Name=stringtie.1.1_m.290
contig_001_length_2062630	transdecoder	five_prime_UTR	9894	9993	.	-	.	ID=stringtie.1.1|m.290.utr5p1;Parent=stringtie.1.1|m.290
contig_001_length_2062630	transdecoder	five_prime_UTR	9745	9746	.	-	.	ID=stringtie.1.1|m.290.utr5p2;Parent=stringtie.1.1|m.290
contig_001_length_2062630	transdecoder	exon	9894	9993	.	-	.	ID=stringtie.1.1|m.290.exon1;Parent=stringtie.1.1|m.290
...
```

The `label` tag then needs to be updated in the `trackList.json` file. Additionally, the UTR color can be changed in the style category as well. By default, it is a dark purple, here changed to a light blue.

```
 "style" : {
    "className" : "feature",
    "color" : "#005824",
    "label" : "name",
    "utrColor" : "#a6bddb"
 },
```

### Coloring by strand ###
For analysing things like anti-sense transcription (as can be done here, since the RNAseq libraries were strand-specific), it is sometimes useful to color transcripts by strand, with forward and reverse as different colors.

The color of features can be changed for [CanvasFeatures](https://jbrowse.org/docs/canvas_features.html) by editing the `trackList.json` file. The `color` tag is a sub-tag of `style`, and here is [given a function](http://gmod.org/wiki/JBrowse_FAQ#How_do_I_customize_feature_colors_.28with_CanvasFeatures.29), which returns `#4658c3` (blue) if the `strand` is 1 (forward strand), and `#46c385` (teal) for all other cases, 0 (no strand given) or -1 (reverse strand). The single line starting with `"color"` is added inside of the `style` tag, and contains a function.

```
 "style" : {
     "className" : "feature",
     "color" : "function(feature) { return feature.get('strand')==1 ?'#4658c3':'#46c385'; }"
 },
```

Like `strand`, features can be colored by score, or [any other parameter from the GFF](http://gmod.org/wiki/JBrowse_FAQ#How_do_I_access_data_about_my_features_in_my_callback_or_plugin).

Some of the histogram parameters can also be changed, to match color schemes. Here, the histogram color is changed to match the darker features in this track.

```
 "histograms" : {
    "color" : "#005824",
 },
```

### Using AUGUSTUS or BRAKER ###
The output of AUGUSTUS is also not a standard GFF format, therefore needs to be modified.

`flatfile-to-json.pl --gff PlacoH13_BRAKER1_augustus_no_comment.gff --trackType CanvasFeatures  --trackLabel AUGUSTUS --out ./ `

### Adding protein matches and domains ###

```
blast2genomegff.py -b hoilungia_vs_hsapiens_blastp_e-3.tab -p blastp -S -g ../tracks/Hhon_BRAKER1_CDS.gff3 -d ~/db/human_uniprot.fasta -x -G > hoilungia_vs_hsapiens_blastp_braker_cds.gff
flatfile-to-json.pl --gff hoilungia_vs_hsapiens_blastp_braker_cds.gff --trackType CanvasFeatures  --trackLabel blastp_v_human --out ./
```

### Adding synteny blocks ###
Synteny blocks spanning multiple genes can be displayed. First, blast the two protein sets (query species against target species). Here, for simplicity, only the first transcript model (called `t1`) is used. This simplifies the downstream processing.

`blastp -query Hhon_BRAKER1_proteins.fasta -db triad_augustus_t1_only.prot.fasta -outfmt 6 -evalue 1e-3 -num_threads 4 > hoilungia_vs_trichoplax_blastp_e-3.tab`

Generate the microsynteny GFF file using the `microsynteny.py` script, by adding the option `-G`.

`microsynteny.py -b hoilungia_vs_trichoplax_blastp_e-3.tab -q Hhon_BRAKER1_genes.gff3 -d Trichoplax_scaffolds_JGI_AUGUSTUS_transcript_only.gff --blast-query-delimiter "." --blast-db-delimiter __ -G > hoilungia_vs_trichoplax_microsynteny_v2.gff`

This is then converted to `CanvasFeatures`, as above. Name fields should appear in jbrowse as `blk-X_to_scaffold_Y`, where X is the block number from the `microsynteny.py`, and `scaffold_Y` is whatever the name is of the target scaffold (e.g. `scaffold_123`, `contig013579`, `chr9`).

`flatfile-to-json.pl --gff hoilungia_vs_trichoplax_microsynteny_v2.gff --trackType CanvasFeatures --trackLabel Triad_microsynteny --out ./`

### Raw read alignment and coverage ###
As [per the instructions](https://jbrowse.org/docs/tutorial_classic.html#next-gen-reads-bam), `BAM` files can be used directly, as long as the index file is there. Otherwise, the index can be regenerated with `samtools index`.

`~/samtools-1.8/samtools index PlacoH13_final_contigs_bowtie2.bam.sorted.bam`

BAM tracks can be directly added to the configuration file `tracks.conf`, here adding both the genomic DNA coverage, and the RNAseq mapping.

```
[tracks.bowtie2dna]
category = NGS
storeClass = JBrowse/Store/SeqFeature/BAM
urlTemplate = PlacoH13_final_contigs_bowtie2.bam.sorted.bam
type = JBrowse/View/Track/Alignments2
key = bowtie2 DNA reads
[tracks.tophat2rnaseq]
category = NGS
storeClass = JBrowse/Store/SeqFeature/BAM
urlTemplate = PlacoH13_hardmasked_tophat2.bam
type = JBrowse/View/Track/Alignments2
key = tophat2 RNAseq reads
```

Depending on coverage, this may not be convenient to display all of the reads. Instead a coverage plot (in *BigWig* format) can be generated from the `.bam` file. This involves conversion of the `.bam` file to `.bed` format (using [bedtools](https://bedtools.readthedocs.io/en/latest/), and then to `.bw`. The `bedGraphToBigWig` program ([downloaded from UCSC here](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)) is used to generate a `BigWig` file. The binary probably needs to be given executable privileges (`chmod +x bedGraphToBigWig`). This program requries two inputs. One is a bed-graph format file of 4 tab-delimited columns, as `scaffold  start  end  coverage`. This is generated by piping the SAM or BAM file into `bedtools`.

`samtools view -b PlacoH13_final_contigs_bowtie2.bam.sorted.bam | bedtools genomecov -ibam stdin -bga > PlacoH13_final_contigs_bowtie2.bam.sorted.cov`

The output looks like this:

```
contig_001_length_2062630	0	1	19
contig_001_length_2062630	1	2	31
contig_001_length_2062630	2	3	33
```

The other contains two columns, one of the scaffold name, the other of scaffold length in bases. This can be generated with many programs, here I use [sizecutter.py](https://bitbucket.org/wrf/sequences/src/master/sizecutter.py), and pipe the output through `sed` to replace the spaces with tabs.

`sizecutter.py -f ../sequences/Hhon_final_contigs_unmasked.fasta | sed s/" "/"\t"/ > Hhon_final_contigs_unmasked.fasta.sizes`

Create the BigWig fie using `bedGraphToBigWig`.

`bedGraphToBigWig PlacoH13_final_contigs_bowtie2.bam.sorted.cov Hhon_final_contigs_unmasked.fasta.sizes PlacoH13_final_contigs_bowtie2.bw`

The above procedure is then done again for the RNAseq reads, using the options `-split` (meaning do not count intervals where reads are split, i.e. RNAseq reads bridging exons) and `-bg`.

`samtools view -b PlacoH13_hardmasked_tophat2.bam | bedtools genomecov -ibam stdin -bg -split > PlacoH13_hardmasked_tophat2.bam.cov`

Then convert to BigWig exactly as above:

`bedGraphToBigWig PlacoH13_hardmasked_tophat2.bam.cov Hhon_final_contigs_unmasked.fasta.sizes PlacoH13_hardmasked_tophat2.bw`

## Some configuration notes for JBrowse ##
It is sometimes useful or convenient to [set up multiple genomes on a single JBrowse instance](http://gmod.org/wiki/JBrowse_FAQ#How_do_I_set_up_multiple_genomes_in_a_single_jbrowse_instance.3F).

Thus, in the normal folder tree of `/var/www/html/jbrowse/`, one could have additional folders inside `jbrowse/` (e.g. `/var/www/html/jbrowse/data/Nvec`) or even symbolic links, say to another drive (e.g. `/mnt/genome_data`.

```
cd /var/www/html/jbrowse/
ln -s /mnt/genome_data/ data
```

Another drive is used to store the data, rather than the directory `/var/www/`. Note that the default Apache2 user (`www-data`) may need to be configured for certain folders, or given `+x` permissions (i.e. of the entire directory tree up to the data directory). In this case, all key files should have `+r` permissions (they probably already do), and directories up to that folder should have `+x`. This would mean that in the above example, the user `www-data` must have permissions to access `/mnt/` and `/mnt/genome-data/`.

