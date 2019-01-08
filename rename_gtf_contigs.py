#!/usr/bin/env python
# v1.0 created 2016-05-16

'''rename_gtf_contigs.py    last modified 2019-01-09

    rename scaffolds/contigs in a GTF/GFF based on a conversion vector

rename_gtf_contigs.py -c conversions.txt -g genes.gtf > renamed_genes.gtf

    conversion vector should be in format of oldgene tab newgene, such as:
oldnamecontig123    newnamecontig001

    or can be in reversed order using -R, as newname tab oldname

    contigs not included in vector are kept as is unless -n is used
    a note is added to the attributes "Rename=false" for removal by grep

    using conversion vector list from:
number_contigs_by_length.py -c conversions.txt contigs.fasta > renamed_contigs.fasta
'''

import sys
import time
import argparse

def make_conversion_dict(conversionfile, do_reverse):
	'''return dict where keys are old contig names and values are new contig names'''
	conversiondict = {}
	print >> sys.stderr, "# Reading conversion file {}".format(conversionfile)
	for line in open(conversionfile,'r'):
		line = line.strip()
		if line:
			if do_reverse:
				conversiondict.update( dict( [(line.split('\t'))[::-1]] ) )
			else:
				conversiondict.update( dict( [(line.split('\t'))] ) )
	print >> sys.stderr, "# Found names for {} contigs".format(len(conversiondict))
	return conversiondict

def make_exclude_dict(excludefile):
	print >> sys.stderr, "# Reading exclusion list {}".format(excludefile), time.asctime()
	exclusion_dict = {}
	for term in open(excludefile,'r'):
		term = term.rstrip()
		if term[0] == ">":
			term = term[1:]
		exclusion_dict[term] = True
	print >> sys.stderr, "# Found {} contigs to exclude".format(len(exclusion_dict) ), time.asctime()
	return exclusion_dict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-c','--conversion', help="text file of naming conversions", required=True)
	parser.add_argument('-E','--exclude', help="file of list of bad contigs")
	parser.add_argument('-g','--gtf', help="gtf or gff format file", required=True)
	parser.add_argument('-n','--nomatch', action="store_true", help="exclude features with no conversion")
	parser.add_argument('-R','--reversed', action="store_true", help="conversion vector is in reversed order, as newname--oldname")
	args = parser.parse_args(argv)

	conversiondict = make_conversion_dict(args.conversion, args.reversed)

	exclusiondict = make_exclude_dict(args.exclude) if args.exclude else None

	linecounter = 0
	conversions = 0
	noconvertfeatures = 0
	print >> sys.stderr, "# Reading features from {}".format(args.gtf), time.asctime()
	for line in open(args.gtf,'r'):
		line = line.strip()
		if line: # remove empty lines
			if line[0]=="#": # write out any comment lines with no change
				print >> wayout, line
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
						print >> sys.stderr, "WARNING: NO CONVERSION FOR {}".format(scaffold)
						lsplits[8] = "{};Rename=false".format(lsplits[8])
				else:
					conversions += 1
					lsplits[0] = newscaffold
				print >> wayout, "\t".join(lsplits)
	print >> sys.stderr, "# Counted {} lines".format(linecounter), time.asctime()
	print >> sys.stderr, "# Converted {} lines and could not change {}".format(conversions, noconvertfeatures), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
