#!/usr/bin/env python

'''reverse_gff_strand.py  last modified 2022-02-11

    reverse all signs from +/- and -/+

reverse_gff_strand.py genes.gff > genes.reversed.gff
'''

import sys
import argparse

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('gff', help="GFF or GTF file")
	parser.add_argument('--types', help="comma-separated list of types of features to reverse, default is all")
	args = parser.parse_args(argv)

	direction_counts = {"+":0, "-":0, ".":0, "?":0}

	accepted_types = []
	if args.types:
		accepted_types = args.types.split(",")
		sys.stderr.write( "# Reversing strands in {}, for features of {}\n".format(args.gff, " ".join(accepted_types) ) )
	else:
		sys.stderr.write( "# Reversing strands in {}\n".format(args.gff) )
	for line in open(args.gff,'r'):
		line = line.strip()
		if line and line[0] != "#":
			lsplits = line.split("\t")
			feature = lsplits[2]
			strand = lsplits[6]
			try:
				direction_counts[strand] = direction_counts.get(strand) + 1
			except KeyError: # print weird line, and add to dict
				sys.stderr.write("WARNING: UNKNOWN STRAND {} on line \n{}\n".format(strand,line) )
				direction_counts[strand] = direction_counts.get(strand,0) + 1
			# if filtering types, and is not the correct type, then print line as is and move on
			if args.types and feature not in accepted_types: # just print as is
				print( line, file=wayout )
				continue

			# if not skipped, then swap strands
			if strand=="+":
				lsplits[6] = "-"
			elif strand=="-":
				lsplits[6] = "+"
			# print updated strand line
			print( "\t".join(lsplits) , file=wayout )
	sys.stderr.write( "# Swapped {} forward and {} reverse, {} with no strand were unchanged\n".format( direction_counts.get("+"), direction_counts.get("-"), direction_counts.get(".") ) )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
