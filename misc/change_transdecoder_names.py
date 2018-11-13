#!/usr/bin/env python

'''convert stringtie GTF to standard GFF (ID-Parent format)
change_transdecoder_names.py transdecoder.genome.gff > transdecoder.genome.renamed.gff

   change transcript features with -t, such as
   -t mRNA
'''

import sys
import argparse

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('gff', help="TransDecoder GFF file")
	args = parser.parse_args(argv)

	geneid_dict = {} # store used IDs, to prevent duplicate features
	print >> sys.stderr, "# Converting Name fields in {}".format(args.gff)
	for line in open(args.gff,'r'):
		line = line.strip()
		if line and line[0] != "#":
			lsplits = line.split("\t")
			feature = lsplits[2]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field])
			if feature=="gene": # changed to ID and Name
				geneid = attrd.get("ID")
				newattrs = "ID={0};Name={0}".format( geneid )
				lsplits[8] = newattrs
				if geneid not in geneid_dict: # should remove duplicate gene IDs
					# in case listed more than once, probably for each mRNA
					geneid_dict[geneid] = True
				else:
					continue
			elif feature=="mRNA": # changed to ID and Name
				try: # for version 3+
					geneid = attrd.get("ID").split("::")[1]
				except IndexError: # for versions 2.x
					geneid = attrd.get("ID").replace("|","_")
				newattrs = "ID={0};Parent={1};Name={2}".format(attrd.get("ID"),attrd.get("Parent"),geneid)
				lsplits[8] = newattrs
			print >> wayout, "\t".join(lsplits)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
