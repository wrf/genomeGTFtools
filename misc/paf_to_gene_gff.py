#!/usr/bin/env python

# paf_to_gene_gff.py 2020-02-26

'''paf_to_gene_gff.py  last modified 2020-02-26

  convert minimap PAF format file to a rough approximation of GFF
  this is meant to get gene coordinates from mapped transcripts
  as was done with the Pleurobrachia bachei genome

  to accurately preserve exon structures, it may be better to use GMAP

paf_to_gene_gff.py pbachei_03_filtered_gene_models_transcripts.paf > pbachei_03_filtered_gene_models_transcripts.rough_bounds.gff

'''

import sys

# PAF is a text format describing the approximate mapping positions between two set of sequences.
#Col 	Type 	Description
#1 	string 	Query sequence name
#2 	int 	Query sequence length
#3 	int 	Query start (0-based; BED-like; closed)
#4 	int 	Query end (0-based; BED-like; open)
#5 	char 	Relative strand: "+" or "-"
#6 	string 	Target sequence name
#7 	int 	Target sequence length
#8 	int 	Target start on original strand (0-based)
#9 	int 	Target end on original strand (0-based)
#10 	int 	Number of residue matches
#11 	int 	Alignment block length
#12 	int 	Mapping quality (0-255; 255 for missing)

#                          query   len start end strand       target     len    start   end   match bllen qual   misc
#scaffold3_1_size219608_gene_7Barcelona	402	9	196	-	AVPN01000021.1	175819	1840	2028	106	188	4	tp:A:P	cm:i:12	s1:i:106	s2:i:93	dv:f:0.0965	rl:i:87
#scaffold3_1_size219608_gene_8Barcelona	5673	34	1981	-	AVPN01001762.1	18701	3731	8325	1538	4594	30	tp:A:P	cm:i:421	s1:i:1477	s2:i:1295	dv:f:0.0209	rl:i:266
#scaffold3_1_size219608_gene_8Barcelona	5673	1075	3027	+	AVPN01000195.1	79835	24426	29571	128	5145	39	tp:A:P	cm:i:25	s1:i:117	s2:i:89	dv:f:0.2132	rl:i:266
#scaffold3_1_size219608_gene_9Barcelona	9642	2848	5045	-	AVPN01000310.1	63752	15433	49119	389	33686	13	tp:A:P	cm:i:126	s1:i:370	s2:i:343	dv:f:0.1161	rl:i:218
#scaffold3_1_size219608_gene_9Barcelona	9642	648	1819	-	AVPN01000026.1	164099	99923	109771	191	9849	60	tp:A:P	cm:i:44	s1:i:177	s2:i:93	dv:f:0.1446	rl:i:218

if len(sys.argv) < 2:
	sys.exit( __doc__ )
else:
	used_txs = {} # key is transcript, value is bool
	sys.stderr.write( "# Reading {}\n".format( sys.argv[1] ) )
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			transcript_name = lsplits[0]
			if used_txs.get(transcript_name,False):
				continue
			used_txs[transcript_name] = True
			scaffold = lsplits[5]
			strand = lsplits[4]
			tstart = lsplits[7]
			tend = lsplits[8]
			match = lsplits[9]
			shortname = transcript_name.split("_")[-1]
			#          scaf                st end strand score  attr
			outline = "{}\tminimap2\tgene\t{}\t{}\t{}\t{}\t.\tID={};Name=gene_{}\n".format( scaffold, tstart, tend, match, strand, transcript_name, shortname )
			sys.stdout.write( outline )
	sys.stderr.write( "# Wrote {} transcripts\n".format( len(used_txs) ) )


#
