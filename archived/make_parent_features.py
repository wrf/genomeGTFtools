#!/usr/bin/env python

"""make_parent_features.py  last modified 2022-04-26
    make parent features for GFF-type files, based on shared ID= tag

    for subfeatues, ID= tag will be renamed to Parent= , omitting the ID= tag

    script will name toplevel feature as 'mRNA', and sublevel as 'exon'

make_parent_features.py features.gff > features_w_parent.gff

    for example:

./make_parent_features.py test_data/Hhon_TRINITY_transcripts_GMAP.gff3 > test_data/Hhon_TRINITY_transcripts_GMAP.w_parent.gff3

"""

import sys
import re
import time
from collections import defaultdict,Counter

if len(sys.argv) < 2:
	sys.stderr.write( __doc__ )
else:
	toplevel = "mRNA"
	sublevel = "exon"

	ID_order = defaultdict(dict) # key is scaffold, then ID, then start position
	feature_dict = defaultdict(list)
	linecounter = 0
	sys.stderr.write("# Reading features from {}  {}\n".format(sys.argv[1], time.asctime() ) )
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line and line[0]!="#":
			linecounter += 1
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			program = lsplits[1]
			lsplits[2] = sublevel # reassign feature type
			start = int(lsplits[3])
			attributes = lsplits[8]
			if attributes.find("ID") >= 0: # gff3 format but no ID
				geneid = re.search("ID=([\w.|-]+)", attributes).group(1)
			lsplits[8] = lsplits[8].replace("ID=","Parent=")
			if geneid not in ID_order[scaffold]:
				ID_order[scaffold][geneid] = start
			else:
				if ID_order[scaffold].get(geneid) > start:
					ID_order[scaffold][geneid] = start
			feature_dict[geneid].append(lsplits)
	sys.stderr.write("# Counted {} lines for {} features\n".format( linecounter, len(feature_dict) ) )

	writecounter = 0
	sys.stderr.write("# Generating parent features  {}\n".format( time.asctime() ) )
	for scaffold in sorted(ID_order.keys()):
		for geneid, pos in sorted(ID_order[scaffold].items(), key=lambda x: x[1]):
			# for a few parameters, take the most extreme of all subfeatures
			# including bounds, score, and also strand
			allstarts = map(int,[ls[3] for ls in feature_dict[geneid]])
			allends = map(int,[ls[4] for ls in feature_dict[geneid]])
			parentstart = min(allstarts)
			parentend = max(allends)
			topscore = max(map(int,[ls[5] for ls in feature_dict[geneid]]))
			# parent feature takes most common strand, which should all be the same
			mainstrand = Counter([ls[6] for ls in feature_dict[geneid]]).most_common(1)[0][0]
			# make and write parent line
			outline = "{0}\t{1}\t{7}\t{2}\t{3}\t{4}\t{5}\t.\tID={6};Name={6}\n".format( scaffold, program, parentstart, parentend, topscore, mainstrand, geneid, toplevel )
			sys.stdout.write( outline )
			writecounter += 1
			# write lines for subfeatures
			for outsplits in feature_dict.get(geneid):
				outline = "{}\n".format( "\t".join(outsplits) )
				sys.stdout.write( outline )
	sys.stderr.write("# Wrote {} parent features\n".format( writecounter ) )
#
