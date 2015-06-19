#!/usr/bin/env python

# v1.0 created 2015-06-16
#
# groupestexons.py use EST information to group exons into plausible transcripts
#
# and using information from:
# https://www.sanger.ac.uk/resources/software/gff/spec.html
# http://www.sequenceontology.org/gff3.shtml

'''
GROUPESTEXONS.PY v1 2015-06-19

  ## GENERAL PRINCIPLES
  assume EST information is always right, or more right than de novo genes
  assume that different ESTs with the same exons are the same transcript
  alternate exon use can be defined by the EST
  ignore antisense transcription
  within a forward reverse pair, any number of exons can be missing

groupestexons.py -e ests.gmap.gff -g gene_models.gff

  specify the delimiter for the ESTs and reads with -n
  default is ., such as for seq123.1 and seq123.2

  ## GETTING GFF RESULTS WITH GMAP

  GMAP should use -f 3 (format 3 gff) as output
  such as:

gmap -d monbr1 Monbr1_cDNA_reads.fasta -f 3 > Monbr1_cdna_reads.gmap.gff

  use multithread option -t for faster results
'''

# TEST PARAMETERS

# TEST CASES
# scaffold_17 167000 XYM20158.x1 forward and reverse overlap, one blast hit as only evidence
# scaffold_17 130200 fgenesh2 170022 fake protein, contains exons going the wrong direction
#    real transcripts are bounded on both ends by ESTs
# scaffold_17 125050 short mapping, appx 20bp to provide false fusion
# scaffold_17 67000 run-on protein containing probably several genes and no RNAseq
# scaffold_17 62000 long EST antisense of several exons, bridging different blast hits
# scaffold_17 43000 two long multiexon genes, one with a single blast hit, other with no evidence

#
import sys
import os
import re
import argparse
import time
from collections import defaultdict,namedtuple
from itertools import groupby

### GFF FORMAT
# GFF3 format is a tabular of:
#
# seqid source type start end score strand phase attributes
#
# in the context of blast results, this is converted as:
# seqid is the subject ID, so sseqid
# source is BLAST
# type could be mRNA, exon, CDS, but in this case HSP as the BLAST result
# start must always be less or equal to end, so if BLAST hits the reverse 
#   complement, these must be switched
# score is the bitscore
# strand is either + or - for forward and reverse hits
# phase is intron phase used only for "CDS" type, otherwise is "."
# attributes is a list of tag=value; pairs usually including query ID

def exon_nums_to_set(exondict):
	for k,v in exondict.iteritems():
		exondict[k] = set(v)

def get_read_bounds(exondict, fnc, ind):
	# returns dictionary of mins for forward reads and maxs for reverse reads
	bounds_d = {}
	for k,v in exondict.iteritems():
		# operates like min( position in tuple for x in interval)
		bounds_d[k] = fnc(x[ind] for x in v)
	return bounds_d

def numbers_to_ranges(numberset):
	# makes groups of value-positions, where values in a row end up as the same group
	for i,l in groupby(enumerate(sorted(numberset)), lambda (x,y): y-x):
		l = list(l)
		# returns generator of first and last positions in the group
		yield l[0][1], l[-1][1]

def make_exon_ranges(exondict):
	exoncounter = 0
	rangeDict = {}
	for k,v in exondict.iteritems():
		exonranges = list(numbers_to_ranges(v))
		rangeDict[k] = exonranges
	# dictionary of lists, where lists contain tuples of exon boundaries
	return rangeDict

def within_exon_range(targetexon, lower, upper):
	return lower <= targetexon <= upper

def get_strand_consensus(strandlist, estscaf):
	strandconsensus = float(sum(strandlist))/len(strandlist)
	if strandconsensus > 0:
		return "+"
	elif strandconsensus < 0:
		return "-"
	else:
		print >> sys.stderr, "WARNING EST {} on {} composed of mixed forward and reverse".format(*estscaf)
		return "."

def get_best_part_exon(refgenedict, exon, verbose, exnum, pos, antipos):
	'''for non exact matches, find closest part of an exon'''
	partmatch = namedtuple("partialmatch", "mtype, gexc, extension")
	pm = partmatch(mtype='n', gexc=None, extension=0)
	for gr in refgenedict:
		opm = partmatch(mtype='n', gexc=None, extension=0)
		if exon[pos]==gr[pos]: # pos is either 0 or 1, for start or end of exon or gr
			opm = partmatch(mtype='p', gexc=gr, extension=exon[pos]-gr[antipos])
			if verbose:
				print >> sys.stderr, "Exon {} partial {} at {} {}".format(exnum, exon, gr, refgenedict[gr] )
		elif exon[0]>gr[0] and exon[1]<gr[1]: # for embedded exon
			opm = partmatch(mtype='m', gexc=gr, extension=exon[pos]-gr[antipos])
			if verbose:
				print >> sys.stderr, "Exon {} embedded {} at {} {}".format(exnum, exon, gr, refgenedict[gr] )
		if opm.gexc: # prioritize matched exons, then partial, then embedded
			if pm.gexc:
				if pm.mtype=='m' and opm.mtype=='p': # take partial over embedded
					pm = partmatch(*opm)
				elif pm.mtype=='p' and opm.mtype=='p': # take longer partials
					if pm.extension < opm.extension:
						pm = partmatch(*opm)
			else: # if no pm is found yet, take the first match
				pm = partmatch(*opm)
	# return tuple of best matching exon, or None
	return pm.gexc

def find_matching_exons(exbyreaddict, read, refgenedict, verbose, doreverse):
	'''for each EST exon find any matching gene exons'''
	# each kept exon is in format of (exon to call line, new attribute ID, possible new coords)
	keptexons = []
	# establish order of first, last exons
	order = (0,len(exbyreaddict[read])-1)
	if doreverse:
		order = tuple(reversed(order))
	# iterate through exons in sorted order of the direction of the sequence (+/-)
	for i,exon in enumerate(sorted(exbyreaddict[read],key=lambda x: x[0], reverse=doreverse) ):
		refex = refgenedict.get(exon, False)
		if refex:
			if verbose:
				print >> sys.stderr, "Exon {} found {} {}".format(i, exon, refex)
			keepexon = (exon, read, exon)
			keptexons.append(keepexon)
		else:
			### TODO add case for good EST but no de novo model
			if verbose:
				print >> sys.stderr, "No exact match for {}".format(exon)
			if i==order[0]: # first exon, 5prime
				bestmatch = get_best_part_exon(refgenedict, exon, verbose, exnum=i, pos=1, antipos=0)
				if bestmatch:
					keepexon = (bestmatch, read, exon) # take positions from EST
					keptexons.append(keepexon)
			# not elif for case of only one exon
			if i==order[1]: # last exon, 3prime
				bestmatch = get_best_part_exon(refgenedict, exon, verbose, exnum=i, pos=0, antipos=1)
				if bestmatch:
					if doreverse:
						keepexon = (bestmatch, read, (bestmatch[0], exon[1]) ) # reverse strand
					else:
						keepexon = (bestmatch, read, (exon[0], bestmatch[1]) ) # forward strand
					keptexons.append(keepexon)
	return keptexons

def print_updated_gff(matchingexons, ref_lines, scaffold, transstrand, program):
	# generate transcript line
	transstart = str(min(me[2][0] for me in matchingexons) )
	transend = str(max(me[2][1] for me in matchingexons) )
	transattrs = 'gene_id "{0}"; transcript_id "{0}.1";'.format(matchingexons[0][1])
	translist = [scaffold, program, "transcript", transstart, transend, "100", transstrand, ".", transattrs]
	print >> sys.stdout, "\t".join(translist)
	exoncount = 0
	# me is a tuple of matching ref exon, read for updating attributes, and possible new positions
	# set is used to remove duplicate exons from forward and reverse
	for me in sorted(set(matchingexons), key=lambda x: x[2][0]):
		exoncount += 1
		outsplits = ref_lines[me[0]].split("\t")
		outsplits[1] = program
		outsplits[3] = str(me[2][0]) # from tuple like (1160, 1498) should be 1160
		outsplits[4] = str(me[2][1]) # and 1498
		attributes = '{} exon_number "{}";'.format(transattrs, exoncount)
		outsplits[8] = attributes
		outline = "\t".join(outsplits)
		print >> sys.stdout, outline

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-e','--ests', help="EST gmap results file")
	parser.add_argument('-g','--genes', help="gene models that will be clustered or corrected")
	parser.add_argument('-b','--blast', help="optional tabular blast for direction evidence")
	parser.add_argument('-l','--map-length', type=int, default=100, help="minimum mapping length in bases [100]")
	parser.add_argument('-i','--intron-distance', type=int, default=1000, help="max distance for introns [1000]")
	parser.add_argument('-r','--read-limit', type=int, default=10, help="max allowed reads per EST group [10]")
	parser.add_argument('-n','--number-split', default=".", help="delimited for EST numbers [.]")
	parser.add_argument('-p','--program', help="program for 2nd column in output [EST]", default="EST")
	parser.add_argument('-v','--verbose', help="verbose output", action="store_true")
	args = parser.parse_args(argv)

	# counter for number of lines, and inferred exons
	linecounter = 0
	exoncounter = 0
	# counter for number of hits that are filtered
	badhits = 0

	starttime = time.time()

	# PART 1
	#
	# make this a dictionary of lists, containing a tuple of exon start and stop positions
	print >> sys.stderr, "Starting exon parsing on %s" % (args.genes), time.asctime()
	gene_ref_exons = defaultdict(dict) # keys are scaffolds, then tuples
	gene_ref_lines = defaultdict(dict) # keys are scaffolds, then lists from lsplits
	for line in open(args.genes,'r'):
		linecounter += 1
		line = line.rstrip()
		if line and not line[0]=="#":
			lsplits = line.split("\t")
			if lsplits[2]=="exon":
				exoncounter += 1
				scaffold = lsplits[0]
				strand = lsplits[6]
				exonbound = ( int(lsplits[3]),int(lsplits[4]) )
				if strand=='+':
					gene_ref_exons[scaffold][exonbound] = 1
				else: # assume '-' strand
					gene_ref_exons[scaffold][exonbound] = -1
				# keep track of all lists for later output
				gene_ref_lines[scaffold][exonbound] = line
			elif lsplits[2]=="CDS":
				pass ### TODO maybe this is needed
	print >> sys.stderr, "Parsed %d lines" % (linecounter), time.asctime()
	print >> sys.stderr, "Counted %d putative exons" % (exoncounter), time.asctime()

	# OPTIONAL
	#
	if args.blast:
		blasthitcount = 0
		blasthitdict = defaultdict(dict)
		for line in open(args.blast, 'r'):
			blasthitcount += 1
			qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.rstrip().split("\t")
			# convert strings of start and end to integers for calculations
			isend = int(send)
			isstart = int(sstart)
			iqend = int(qend)
			iqstart = int(qstart)

			# is checks are passed, accumulate that hit an info into dict
			hittuple = (iqstart, iqend, isstart, isend)
			# must use setdefault since get(qseqid, []) returns None, and cannot append
			blasthitdict[sseqid].setdefault(qseqid, []).append(hittuple)

	# PART II
	#
	# read EST mappings
	plusexons, minusexons = 0,0	
	plusstrand, minusstrand = 0,0
	# counter for potentially overlapping hits or tandem duplicates
	tandemdups = 0
	hitDictCounter = defaultdict(int)

	#print >> sys.stderr, "Reading database sequences from %s" % (args.db), time.asctime()
	#dbdict = SeqIO.to_dict(SeqIO.parse(args.db, 'fasta'))
	#print >> sys.stderr, "Counted %d db sequences" % (len(dbdict)), time.asctime()

	#est_exons = defaultdict(list)
	estexoncount = 0
	readlengths = defaultdict(int) # sum of exons by readname
	estparts = defaultdict(list) # list of two possible readnames for each estname
	exonsbyscaffold = defaultdict( lambda: defaultdict(list) ) # by scaffold then readname
	estbyread = {}
	forwardreads = defaultdict(list)
	reversereads = defaultdict(list)
	#scafreadcount = defaultdict( lambda: defaultdict(int) )

	readidre = "Name=([\w\d\.|:]+);" # must have . | and : for different EST types
	print >> sys.stderr, "Reading EST sequences from %s" % (args.ests), time.asctime()
	for line in open(args.ests, 'r'):
		estexoncount += 1
		line = line.rstrip()
		if line and not line[0]=="#":
			lsplits = line.split("\t")
			attrs = lsplits[8]
			# from ID=XYM14183.y1.path1;Name=XYM14183.y1;Target=XYM14183.y1 88 130;Gap=M43
			# or ID=jgi|JGI_XYM10825.rev|.path1;Name=jgi|JGI_XYM10825.rev|;Target=jgi|JGI_XYM10825.rev|
			# or ID=3720288:1.path1;Name=3720288:1;Target=3720288:1 12 684;Gap=M673
			readname = re.search(readidre,attrs).group(1) # should be XYM14183.y1
			estname = readname.split(args.number_split)[0] # should be XYM14183

			### TODO parse this better to allow for multiple paths as different genes
			estbyread[readname] = estname
			estparts[estname].append(readname)
			scaffold = lsplits[0]
			strand = lsplits[6]
			exonbound = ( int(lsplits[3]),int(lsplits[4]) )

			exonsbyscaffold[scaffold][readname].append(exonbound)
			#est_exons[scaffold].extend(range(*exonbound) ) # ests are treated strand-non-specific
			readlengths[readname] += exonbound[1]-exonbound[0]+1
			#scafreadcount[scaffold][readname] += 1
			if strand=='+':
				forwardreads[readname].append(exonbound)
			else: # assume '-' strand
				reversereads[readname].append(exonbound)
	print >> sys.stderr, "Counted %d forward and %d reverse reads" % (len(forwardreads), len(reversereads) ), time.asctime()
	print >> sys.stderr, "Counted %d ESTs with %d bases" % (len(readlengths), sum(readlengths.values()) ), time.asctime()
	print >> sys.stderr, "Counted %d putative exons" % (estexoncount), time.asctime()

	forwardbounds = get_read_bounds(forwardreads, fnc=min, ind=0)
	reversebounds = get_read_bounds(reversereads, fnc=max, ind=1)

	# if any exons are present from one read in another, take the exons and pop the read

	pairscount = 0
	for k,v in estparts.iteritems():
		if all(bool(r) for r in v):
			pairscount+=1
	print >> sys.stderr, "Counted %d ESTs with aligned pairs" % (pairscount), time.asctime()

	#print >> sys.stderr, "Merging ESTs", time.asctime()
	#exon_nums_to_set(est_exons)
	#merged_exons = make_exon_ranges(est_exons)
	#print >> sys.stderr, merged_exons
	#print >> sys.stderr, "Extending ranges of ESTs", time.asctime()

	matchcount = 0
	transcriptcount = 0
	print >> sys.stderr, "Finding predicted exons that match ESTs", time.asctime()
	for sc,ebr in exonsbyscaffold.iteritems(): # scaffold, dicts of exons by read
		#extendedranges = {}
		#readsbyexon = defaultdict(set)
		#exbyread = defaultdict(list)
		#for read, exs in ebr.iteritems(): # read, exons
		#	for ex in exs:
		#		for me in merged_exons[sc]:
		#			if within_exon_range(ex[0], *me) or within_exon_range(ex[1], *me):
		#				readsbyexon[me].add(read)
		#				exbyread[read].append(me)

		readgroups = defaultdict(list) # should be list by final transcript boundaries
		readonscaffold = ebr.keys()
		estonscaffold = list(set(estbyread[r] for r in readonscaffold))
		for est in estonscaffold:
			forwardmatches = []
			reversematches = []
			matchingexons = []
			if args.verbose:
				print >> sys.stderr, "Sorting {} on scaffold {}".format(est, sc)
		#	if not est=="3711282":
		#		continue
			readsperest = list(set(estparts[est]))
			if len(readsperest) > args.read_limit:
				print >> sys.stderr, "WARNING counted {} reads from {}, skipping".format(len(readsperest), est)
				continue
			for read in readsperest:
				if args.verbose:
					print >> sys.stderr, "Sorting exons from {}".format(read)
				# matches are in format of (exon, read, new values for position)
				forwardmatches.extend(find_matching_exons(forwardreads, read, gene_ref_exons[sc], args.verbose, doreverse=False) )
				reversematches.extend(find_matching_exons(reversereads, read, gene_ref_exons[sc], args.verbose, doreverse=True) )
			try:
				forward_end = max(me[2][1] for me in forwardmatches)
			except ValueError: # for empty list max fails, so no forward read
				forward_end = 0
			try:
				reverse_end = min(me[2][0] for me in reversematches)
			except ValueError: # for empty list min fails, meaning no reverse read
				reverse_end = 0
			matchingexons.extend(forwardmatches)
			matchingexons.extend(reversematches)

			# add exons in between two ends of reads
			### TODO add optional boundaries for cases with only one read sense
			if forward_end and reverse_end:
				if reverse_end > forward_end: # only if reverse_end is greater than forward_end
					if args.verbose:
						print >> sys.stderr, "Checking internal exons between {} and {}".format(forward_end, reverse_end)
					betweenmatches = []
					for exon in gene_ref_exons[sc]:
						if exon[0]>forward_end and exon[1]<reverse_end:
							### TODO also check against intron distance
							if args.verbose:
								print >> sys.stderr, "Internal exon found {}".format(exon)
							midmatch = (exon, est, exon)
							betweenmatches.append(midmatch)
					matchingexons.extend(betweenmatches)

			if matchingexons:
				transcriptcount += 1
				matchcount += len(matchingexons)
				# forces '+' if there is no gene model
				### TODO incorporate evidence from blast here
				transstrand = get_strand_consensus([gene_ref_exons[sc].get(x[0],1) for x in matchingexons ], (est,sc) )
				print_updated_gff(matchingexons, gene_ref_lines[sc], sc, transstrand, args.program)
	print >> sys.stderr, "Counted %d matching exons or partial exons" % (matchcount), time.asctime()
	print >> sys.stderr, "Printed %d transcripts" % (transcriptcount), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
