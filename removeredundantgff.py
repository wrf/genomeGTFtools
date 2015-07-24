#!/usr/bin/env python
# removeredundantgff.py
# v1.0 2015-07-23

'''
removeredundantgff.py last modified 2015-07-24
    remove identical gene predictions from gff3 file

    EXAMPLE USAGE:

removeredundantgff.py -g genewise.gff3 > genewise.nonredundant.gff3

    IF EXONS WERE PREVIOUSLY EXCLUDED, GET RANGE INFORMATION FROM CDS:

removeredundantgff.py -g genewise.gff3 -C > genewise.nonredundant.gff3
'''

import sys
import os
import argparse
import time
import re
from collections import defaultdict
from itertools import chain,groupby
from numpy import histogram,arange

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-g','--gff', help="gff3 format file")
	parser.add_argument('-d','--delimiter', default=".", help="delimiter for separating names and IDs [.]")
	parser.add_argument('-C','--cds', action="store_true", help="ranges are by CDS rather than exons")
	parser.add_argument('-E','--exons', action="store_false", help="exclude exons when writing output")
	parser.add_argument('-v','--verbose', help="verbose output", action="store_true")
	args = parser.parse_args(argv)

	exonScaffold = defaultdict(lambda: defaultdict(list)) # dict of dicts of lists, as scaffold[gene][exons]
	geneorder = [] # keep order of genes read in the file, to make output deterministic
	scaffoldbygene = defaultdict(list) # all scaffold hits by each gene, might only be one each
	linesbyid = defaultdict(list) # this stores all GFF output lines by the gene ID
	genecount, mrnacount, exoncount, cdscount = 0,0,0,0
	print >> sys.stderr, "Starting exon parsing on {}".format(args.gff), time.asctime()
	for line in open(args.gff,'r').readlines():
		line = line.rstrip()
		if not line or line[0]=="#":
			continue
		lsplits = line.split("\t")
		if lsplits[2]=="gene":
			genecount += 1
			transid = re.search("ID=([\w.]+);", lsplits[8]).group(1)
			geneorder.append(transid)
			scaffoldbygene[transid].append( lsplits[0] )
			linesbyid[transid].append(line)
		elif lsplits[2]=="mRNA":
			mrnacount += 1
			linesbyid[transid].append(line)
		elif lsplits[2]=="exon":
			exoncount += 1
			if not args.cds:
				contig, startpos, endpos, strand = lsplits[0], int(lsplits[3]), int(lsplits[4]), lsplits[6]
				boundaries = (startpos,endpos)
				#if strand=="+":
				exonScaffold[contig][transid].append(boundaries)
				#else: # strand=="-":
			if args.exons:
				linesbyid[transid].append(line)
		elif lsplits[2]=="CDS":
			cdscount += 1
			if args.cds:
				contig, startpos, endpos, strand = lsplits[0], int(lsplits[3]), int(lsplits[4]), lsplits[6]
				boundaries = (startpos,endpos)
				exonScaffold[contig][transid].append(boundaries)
			linesbyid[transid].append(line)
	print >> sys.stderr, "Counted {} gene and {} mRNA predictions".format(genecount, mrnacount), time.asctime()
	print >> sys.stderr, "Counted {} exons and {} CDS".format(exoncount, cdscount), time.asctime()

	printcount = 0
	exonsorted = defaultdict(dict) # dict of scaffold and then exon boundaries as list of tuples
	print >> sys.stderr, "Starting duplicate gene sorting", time.asctime()
	for gene in geneorder:
		for sc in scaffoldbygene[gene]:
			# this generates a string of the sorted list of exon boundaries
			exonorder = str(sorted(exonScaffold[sc][gene]))
			# such as "[(19899, 20013), (20080, 20137), (20218, 20335)]"
			### TODO allow removal of substring matches
			if not exonorder in exonsorted[sc]:
				printcount += 1
				exonsorted[sc][exonorder] = True
				for line in linesbyid[gene]:
					print >> wayout, line
	print >> sys.stderr, "Printed {} non redundant predictions".format(printcount), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
