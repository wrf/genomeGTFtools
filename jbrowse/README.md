# genomeGTFtools

## Tips for use in [JBrowse](https://jbrowse.org/) ##
JBrowse is a genome browser that can run in a web browser. Because of the lack of consistency in GTF/GFF formats (all of which hinges on the flexibility of the final column), various genome browsers will interpret GFF files differently.

Thus, to best configure GFF files for display in JBrowse, there are several considerations.

* Parent features need to be generated for anything in the GFF spanning more than 1 line. This is for memory considerations in how the information is fetched. This means that several of the scripts will have to generate additional lines for a parent feature, rather than relying on the ID of each feature.
* Some features are given by particular programs, but relationships between them are not recognized by default by JBrowse. For instance, [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) will output a GFF with the features `gene-transcript-intron-CDS`. This needs to be converted to `gene-mRNA-CDS` in order to display correctly.
* Older GTF versions will likely need substantial format conversion to display as expected.

## An example workflow ##
Here I make use of some examples from the [genome](https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/) of the [placozoan *Hoilungia hongkongensis*](https://doi.org/10.1371/journal.pbio.2005359).

In preparing files for JBrowse, the default directory for the objects is `/var/www/html/jbrowse/data/`. Here, to set up subfolders to allow multiple genomes in the same JBrowse instance, the local folder is specified in the script `prepare-refseqs.pl` using the option `--out`. 

Another drive is used to store the data, rather than the directory `/var/www/`. Note that the default Apache2 user (`www-data`) may need to be configured for certain folders, or given `+x` permissions (i.e. of the entire directory tree up to the data directory). In this case, all key files should have `+r` permissions (they probably already do), and directories up to that folder should have `+x`.

```
cd /mnt/genome_data/Hhon/
prepare-refseqs.pl --fasta Hhon_final_contigs_unmasked.fasta --out ./
```

Because of issues with parsing the GTF format of [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml), this needs to be converted to a GFF file with a Parent-Child feature, here simply as mRNA-exon. By default, each transcript and exon will be a standalone feature, rather than a transcript composed of exons as subfeatures.

```
stringtie_gtf_to_gff3.py Hhon_tophat2_stringtie.gtf > Hhon_tophat2_stringtie.gff3
flatfile-to-json.pl --gff Hhon_tophat2_stringtie.gff --trackType CanvasFeatures --trackLabel Stringtie --out ./
```

### Coloring by strand ###
For analysing things like anti-sense transcription (as can be done here, since the RNANseq libraries were strand-specific), it is sometimes useful to color transcripts by strand, with forward and reverse as different colors.

The color of features can be changed for [CanvasFeatures](https://jbrowse.org/docs/canvas_features.html) by editing the `trackList.json` file. The `color` tag is a sub-tag of `style`, and here is [given a function](http://gmod.org/wiki/JBrowse_FAQ#How_do_I_customize_feature_colors_.28with_CanvasFeatures.29), which returns `#4658c3` (blue) if the `strand` is 1 (forward strand), and `#46c385` (teal) for all other cases, 0 (no strand given) or -1 (reverse strand). The single line starting with `"color"` is added inside of the `style` tag, and contains a function.

```
"style" : {
    "className" : "feature",
    "color" : "function(feature) { return feature.get('strand')==1 ?'#4658c3':'#46c385'; }"
},
```

Like `strand`, features can be colored by score, or [any other parameter from the GFF](http://gmod.org/wiki/JBrowse_FAQ#How_do_I_access_data_about_my_features_in_my_callback_or_plugin).

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

## Some configuration notes for JBrowse ##
It is sometimes useful or convenient to [set up multiple genomes on a single JBrowse instance](http://gmod.org/wiki/JBrowse_FAQ#How_do_I_set_up_multiple_genomes_in_a_single_jbrowse_instance.3F).

Thus, in the normal folder tree of `/var/www/html/jbrowse/`, one could have additional folders inside `jbrowse/` (e.g. `/var/www/html/jbrowse/data/Nvec`) or even symbolic links, say to another drive (e.g. `/mnt/genome_data`.

```
cd /var/www/html/jbrowse/
ln -s /mnt/genome_data/ data
```

It appears that apache2 needs permissions to access all folders upstream of a target directory to work. This would mean that in the above example, the apache2 user (`www-data`) has permissions to access `/mnt/` and `/mnt/genome-data/`.
