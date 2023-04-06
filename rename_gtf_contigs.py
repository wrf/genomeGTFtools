#!/usr/bin/env python
# v1.0 created 2016-05-16
# python3 update 2019-09-26
# v1.2 allow gzip 2023-01-25
# v1.3 allow skipping of comment lines 2023-04-05

'''rename_gtf_contigs.py    last modified 2023-04-05

    rename scaffolds/contigs in a GTF/GFF based on a conversion vector

rename_gtf_contigs.py -c conversions.txt -g genes.gtf > renamed_genes.gtf

    conversion vector should be in format of oldgene tab newgene, such as:
oldnamecontig123    newnamecontig001

    or can be in reversed order using -R, as newname tab oldname

    contigs not included in vector are kept as is unless -n is used
    a note is added to the attributes "Rename=false" for removal by grep
    this assumes GFF format - NOT GTF

    using conversion vector list from:
number_contigs_by_length.py -c conversions.txt contigs.fasta > renamed_contigs.fasta
'''

import sys
import os
import time
import gzip
import argparse

def make_conversion_dict(conversionfile, do_reverse):
	'''return dict where keys are old contig names and values are new contig names'''
	conversiondict = {}
	sys.stderr.write("# Reading conversion file {}\n".format(conversionfile) )
	for line in open(conversionfile,'r'):
		line = line.strip()
		if line:
			if do_reverse:
				conversiondict.update( dict( [(line.split('\t'))[::-1]] ) )
			else:
				conversiondict.update( dict( [(line.split('\t'))] ) )
	sys.stderr.write("# Found names for {} contigs\n".format(len(conversiondict)) )
	return conversiondict

def make_exclude_dict(excludefile):
	sys.stderr.write("# Reading exclusion list {}  {}\n".format(excludefile, time.asctime() ) )
	exclusion_dict = {}
	for term in open(excludefile,'r'):
		term = term.rstrip()
		if term[0] == ">":
			term = term[1:]
		exclusion_dict[term] = True
	sys.stderr.write("# Found {} contigs to exclude  {}\n".format(len(exclusion_dict) , time.asctime() ) )
	return exclusion_dict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-c','--conversion', help="text file of naming conversions", required=True)
	parser.add_argument('-E','--exclude', help="file of list of bad contigs")
	parser.add_argument('-g','--gtf', help="input file in gtf or gff format, can be .gz", required=True)
	parser.add_argument('-n','--nomatch', action="store_true", help="exclude features with no conversion")
	parser.add_argument('--omit-comments', action="store_true", help="remove any commented lines #")
	parser.add_argument('-R','--reversed', action="store_true", help="conversion vector is in reversed order, as newname--oldname")
	args = parser.parse_args(argv)

	conversiondict = make_conversion_dict(args.conversion, args.reversed)

	exclusiondict = make_exclude_dict(args.exclude) if args.exclude else None

	linecounter = 0
	conversions = 0
	noconvertfeatures = 0

	if args.nomatch:
		sys.stderr.write("# REMOVING scaffolds without matches: {}\n".format(args.nomatch) )
	if args.omit_comments:
		sys.stderr.write("# REMOVING all comment lines: {}\n".format(args.omit_comments) )

	if args.gtf.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing features from {} as gzipped  {}\n".format(args.gtf, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing features from {}  {}\n".format(args.gtf, time.asctime() ) )
	for line in opentype(args.gtf,'rt'):
		line = line.strip()
		if line: # remove empty lines
			if line[0]=="#": # write out any comment lines with no change
				if args.omit_comments:
					continue
				else:
					print( line, file=wayout )
			else:
				linecounter += 1
				lsplits = line.split("\t")
				scaffold = lsplits[0]
				if exclusiondict and exclusiondict.get(scaffold, False):
					continue # skip everything from this contig anyway
				newscaffold = conversiondict.get(scaffold, None)
				if newscaffold is None:
					noconvertfeatures += 1
					if args.nomatch: # no match, so skip
						continue
					else:
						sys.stderr.write("WARNING: NO CONVERSION FOR {}\n".format(scaffold) )
						lsplits[8] = "{};Rename=false".format(lsplits[8])
				else:
					conversions += 1
					lsplits[0] = newscaffold
				print( "\t".join(lsplits), file=wayout )
	sys.stderr.write("# Counted {} lines  ".format(linecounter) + time.asctime() + os.linesep)
	sys.stderr.write("# Converted {} lines and could not change {}\n".format(conversions, noconvertfeatures) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
