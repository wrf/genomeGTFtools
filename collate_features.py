#!/usr/bin/env python

'''
collate_features.py parent_features.gff child_features.gff > combined.gff

'''

import sys
import re
from collections import defaultdict

if len(sys.argv) < 2:
	sys.stderr.write(__doc__)
else:
	subfeature_dict = defaultdict(list)
	for line in open(sys.argv[2],'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			attributes = lsplits[8]
			if attributes.find("Parent")==0: # gff3 format but no ID
				geneid = re.search('Parent=([\w.|-]+)', attributes).group(1)
				subfeature_dict[geneid].append(line)
			else:
				sys.stderr.write("WARNING CANNOT FIND PARENT ID FOR {}\n".format(attributes) )

	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			attributes = lsplits[8]
			if attributes.find("ID")==0: # gff3 format but no ID
				geneid = re.search('ID=([\w.|-]+)', attributes).group(1)
				sys.stdout.write( line + os.linesep )
				for outline in subfeature_dict.get(geneid):
					sys.stdout.write( outline + os.linesep )
			else:
				sys.stderr.write("WARNING CANNOT FIND PARENT ID FOR {}\n".format(attributes) )

