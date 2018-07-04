#!/usr/bin/env python
# microsynteny.py
# v1.0 2015-10-09

'''microsynteny.py v1.2 last modified 2018-07-04

microsynteny.py -q query.gtf -d ref_species.gtf -b query_vs_ref_blast.tab -E ../bad_contigs -g -D '_' --blast-query-delimiter '.' > query_vs_ref_microsynteny.tab

    FOR GENERATION OF BLAST DATA TO LINK SETS
blastx -query query.fasta -db ref_prots.fasta -outfmt 6 -evalue 1e-5 > query_vs_ref_blast.tab

    large blast.tab files can be gzipped as .gz

    IN CASE OF NO OUTPUT, CHECK THAT -Q AND -D ARE SET CORRECTLY
    from gtf files, gene_id is extracted
    the blast query names must match the gene_id
    for example, if gene_id is avic.12345 and blast query is avic.12345.1
    set --blast-query-delimiter '.'

    OUTPUT is tab delimited text file of 11 columns, consisting of:
query-scaffold  ref-scaffold  block-number
  query-gene  start  end  strand  ref-gene  start  end  strand

    THIS CANNOT DETECT ERRONEOUS FUSION OR SPLITTING OF GENES
    i.e. three collinear genes in the query that are erroneously fused
    in the ref species will still count as a block of three

    to determine the correct minimum -m, check the same dataset
    as randomized queries, using -R
    blocks of 2 occur often, but 3 is rare, so -m 3 is usually sufficient
'''

import sys
import re
import time
import argparse
import random
import gzip
from collections import namedtuple,defaultdict

querygene = namedtuple("querygene", "start end strand")
refgene = namedtuple("refgene", "scaffold start end strand")

def parse_gtf(gtffile, exonstogenes, excludedict, isref=False):
	'''from a gtf, return a dict of dicts where keys are scaffold names, then gene names, and values are gene info as a tuple of start end and strand direction'''
	print >> sys.stderr, "# Parsing {}".format(gtffile), time.asctime()
	if isref: # meaning is db/subject, thus get normal dictionary
		genesbyscaffold = {}
	else:
		genesbyscaffold = defaultdict(dict) # scaffolds as key, then gene name, then genemapping tuple
	if exonstogenes:
		nametoscaffold = {} # in order to get transcript boundaries, store names to scaffolds
		nametostrand = {} # store strand by gene ID
		exonboundaries = defaultdict(list) # make list of tuples of exons by transcript, to determine genes
	for line in open(gtffile).readlines():
		line = line.strip()
		if line and not line[0]=="#": # ignore empty lines and comments
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			if excludedict and excludedict.get(scaffold, False):
				continue # skip anything that hits to excludable scaffolds
			feature = lsplits[2]
			attributes = lsplits[8]
			if feature=="gene" or feature=="transcript" or feature=="mRNA":
				try:
					geneid = re.search('gene_id "([\w.|-]+)";', attributes).group(1)
				except AttributeError: # in case re fails and group does not exist
					geneid = re.search('ID=([\w.|-]+);', attributes).group(1)
				if isref:
					refbounds = refgene(scaffold=lsplits[0], start=int(lsplits[3]), end=int(lsplits[4]), strand=lsplits[6] )
					genesbyscaffold[geneid] = refbounds
				else:
					boundstrand = querygene(start=int(lsplits[3]), end=int(lsplits[4]), strand=lsplits[6] )
					genesbyscaffold[scaffold][geneid] = boundstrand
			elif exonstogenes and feature=="exon":
				try:
					geneid = re.search('gene_id "([\w.|-]+)";', attributes).group(1)
				except AttributeError: # in case re fails and group does not exist
					geneid = re.search('name "([\w.|-]+)";', attributes).group(1)
				nametoscaffold[geneid] = scaffold
				nametostrand[geneid] = lsplits[6]
				exonbounds = (int(lsplits[3]), int(lsplits[4]))
				exonboundaries[geneid].append(exonbounds) # for calculating gene boundaries
	if len(genesbyscaffold) > 0: # even if no-genes was set, this should be more than 0 if genes were in one gtf
		print >> sys.stderr, "# Found {} genes".format(sum(len(x) for x in genesbyscaffold.values()) ), time.asctime()
		return genesbyscaffold
	else: # generate gene boundaries by scaffold
		print >> sys.stderr, "# Estimated {} genes from {} exons".format(len(exonboundaries), sum(len(x) for x in exonboundaries.values() ) ), time.asctime()
		for gene,exons in exonboundaries.iteritems():
			if isref: # make different tuple for reference genes, since they are indexed by name, not scaffold
				refbounds = refgene(scaffold=nametoscaffold[gene], start=min(x[0] for x in exons), end=max(x[1] for x in exons), strand=nametostrand[gene] )
				genesbyscaffold[gene] = refbounds
			else:
				boundstrand = querygene(start=min(x[0] for x in exons), end=max(x[1] for x in exons), strand=nametostrand[gene] )
				genesbyscaffold[nametoscaffold[gene]][gene] = boundstrand
		print >> sys.stderr, "# Found {} genes".format(len(genesbyscaffold) ), time.asctime() # uses len here
		return genesbyscaffold

def parse_tabular_blast(blasttabfile, evaluecutoff, querydelimiter, refdelimiter, maxhits=10):
	'''read tabular blast file, return a dict where key is query ID and value is subject ID'''
	if blasttabfile.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# Parsing tabular blastx output {} as gzipped".format(blasttabfile), time.asctime()
	else: # otherwise assume normal open for fasta format
		opentype = open
		print >> sys.stderr, "# Parsing tabular blastx output {}".format(blasttabfile), time.asctime()
	query_to_sub_dict = defaultdict(dict)
	evalueRemovals = 0
	for line in opentype(blasttabfile, 'r').readlines():
		line = line.strip()
		lsplits = line.split("\t")
		# qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
		queryseq = lsplits[0].rsplit(querydelimiter,1)[0]
		if float(lsplits[10]) > evaluecutoff:
			evalueRemovals += 1
			continue
		if len(query_to_sub_dict[queryseq]) >= maxhits: # too many hits already, skip
			continue
		subjectid = lsplits[1].rsplit(refdelimiter,1)[0]
		query_to_sub_dict[queryseq][subjectid] = True
	print >> sys.stderr, "# Found blast hits for {} query sequences".format( len(query_to_sub_dict) ), time.asctime()
	print >> sys.stderr, "# Removed {} hits by evalue".format(evalueRemovals), time.asctime()
	print >> sys.stderr, "# Names parsed as {} from {}, and {} from {}".format( queryseq,lsplits[0], subjectid,lsplits[1] )
	return query_to_sub_dict

def randomize_genes(refdict):
	'''take the gtf dict and randomize the gene names for all genes, return a similar dict of dicts'''
	genepositions = {} # store gene positions as tuples
	randomgenelist = []
	print >> sys.stderr, "# Randomizing query gene positions", time.asctime()
	for scaffold, genedict in refdict.iteritems(): # iterate first to get list of all genes
		for genename in genedict.keys():
			randomgenelist.append(genename)
			genepositions[genename] = genedict[genename]
	# randomize the list
	random.shuffle(randomgenelist)
	# reiterate in same order, but store random gene names at the same position
	genecounter = 0
	randomgenesbyscaf = defaultdict(dict) # scaffolds as key, then gene name, then genemapping tuple
	for scaffold, genedict in refdict.iteritems(): # iterate again to reassign genes to each scaffold
		for genename, bounds in genedict.iteritems():
			randomgenesbyscaf[scaffold][randomgenelist[genecounter]] = genepositions[genename]
			genecounter += 1
	print >> sys.stderr, "# Randomized {} genes".format(genecounter), time.asctime()
	return randomgenesbyscaf

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast', help="tabular blast output", required=True)
	parser.add_argument('-q','--query-gtf', help="gtf file of query genes", required=True)
	parser.add_argument('-d','--db-gtf', help="gtf file of reference genes", required=True)
	parser.add_argument('-Q','--query-delimiter', help="gene transcript separator for query [.]", default='.')
	parser.add_argument('-D','--db-delimiter', help="gene transcript separator for db [.]", default='.')
	parser.add_argument('--blast-query-delimiter', help="gene transcript separator for blast query [|]", default='|')
	parser.add_argument('--blast-db-delimiter', help="gene transcript separator for blast ref [|]", default='|')
	parser.add_argument('-e','--evalue', type=float, default=1e-4, help="evalue cutoff for post blast filtering [1e-4]")
	parser.add_argument('-E','--exclude', help="file of list of bad contigs")
	parser.add_argument('-g','--no-genes', action="store_true", help="genes are not defined, get gene ID for each exon")
	parser.add_argument('-m','--minimum', type=int, default=3, help="minimum syntenic genes to keep block, must be >=2 [3]")
	parser.add_argument('-s','--span', type=int, default=5, help="max number of skippable genes [5]")
	parser.add_argument('-z','--distance', type=int, default=30000, help="max distance before next gene [30000]")
	parser.add_argument('-R','--randomize', help="randomize positions of query GTF", action="store_true")
	parser.add_argument('-v','--verbose', help="verbose output", action="store_true")
	args = parser.parse_args(argv)

	if args.minimum < 2:
		print >> sys.stderr, "WARNING: MINIMUM COLINEARITY -m MUST BE GREATER THAN 1, {} GIVEN".format(args.minimum)
		print >> sys.stderr, "SETTING MINIMUM COLINEARITY TO 2"
		args.minimum = 2

	if args.exclude:
		print >> sys.stderr, "# Reading exclusion list {}".format(args.exclude), time.asctime()
		exclusionDict = {}
		for term in open(args.exclude,'r').readlines():
			term = term.strip()
			if term[0] == ">":
				term = term[1:]
			exclusionDict[term] = True
		print >> sys.stderr, "# Found {} contigs to exclude".format(len(exclusionDict) ), time.asctime()
	else:
		exclusionDict = None

	### SETUP DICTIONARIES ###
	querydict = parse_gtf(args.query_gtf, args.no_genes, exclusionDict)
	refdict = parse_gtf(args.db_gtf, args.no_genes, exclusionDict, isref=True)
	blastdict = parse_tabular_blast(args.blast, args.evalue, args.blast_query_delimiter, args.blast_db_delimiter)

	### IF DOING RANDOMIZATION ###
	if args.randomize:
		querydict = randomize_genes(querydict)

	### START SYNTENY WALKING ###
	syntenylist = []
	blocknum = 1
	blocklengths = defaultdict(int) # dictionary to keep track of number of blocks of length N
	splitgenes = 0 # counter for number of genes where next gene hits same reference, so query might be split
	basetotal = 0 # counter for block length in bases
	scaffoldgenecounts = defaultdict(int)
	lastmatch = "" # default is empty string
	print >> sys.stderr, "# searching for colinear blocks of at least {} genes, with up to {} intervening genes".format( args.minimum, args.span )
	for scaffold, transdict in querydict.iteritems():
		orderedtranslist = sorted(transdict.items(), key=lambda x: x[1].start) # sort by start position
		genesonscaff = len(orderedtranslist) # keeping track of number of genes by scaffold, for scale
		scaffoldgenecounts[genesonscaff] += 1
		if args.verbose:
			print >> sys.stderr, "# Scanning scaffold {0} with {1} genes".format(scaffold, genesonscaff)
		if genesonscaff < args.minimum: # not enough genes, thus no synteny would be found
			if args.verbose:
				print >> sys.stderr, "# Only {1} genes on scaffold {0}, skipping scaffold".format(scaffold, genesonscaff)
			continue
		accountedgenes = [] # list of genes already in synteny blocks on this contig
		for i,transtuple in enumerate(orderedtranslist):
			walksteps = args.span # genes until drop, walk starts new for each transcript
			# starting from each transcript
			startingtrans = transtuple[0]
			refmatch_dict = blastdict.get(startingtrans,None) # get the blast match
			if refmatch_dict==None: # if no blast match, then skip to next gene
				if args.verbose:
					print >> sys.stderr, "# No blast matches for {}, skipping walk".format(startingtrans)
				continue
			for refmatch in refmatch_dict.iterkeys():
				if refmatch==lastmatch: # unique test to see if blast hit matches previous
					splitgenes += 1
				if startingtrans in accountedgenes: # if gene is already in a synteny block, ignore it
					continue
				syntenylist = [(startingtrans,refmatch)]
				if i < genesonscaff - 1: # this allows for 2 genes left
					if args.verbose:
						print >> sys.stderr, "# Starting walk from gene {}".format(startingtrans)
					refcontig = refdict[refmatch].scaffold
					endpos = transtuple[1].end
					# begin of gene walk
					for ntrans, ngeneinfo in orderedtranslist[i+1:]: # getting next transcript, and next gene
						if walksteps > 0:
							if ngeneinfo.start-endpos > args.distance: # next gene is too far
								if args.verbose:
									print >> sys.stderr, "# Next gene {} bases away from {}, stopping walk".format(ngeneinfo.start-endpos, ntrans)
								break # end gene block
							nmatch_dict = blastdict.get(ntrans,None)
							if nmatch_dict==None:
								if args.verbose:
									print >> sys.stderr, "# No blast matches for {}, skipping gene".format(ntrans)
								walksteps -= 1 # if no blast match, then skip to next walk step, and decrement
								continue
							for nmatch in nmatch_dict.iterkeys():
								if nmatch==refmatch: # if last match is the same as current match
									accountedgenes.append(ntrans)
									continue # since it is still the same gene, move on but do not penalize
								elif nmatch.rsplit(args.db_delimiter,1)[0]==refmatch.rsplit(args.db_delimiter,1)[0]:
									accountedgenes.append(ntrans)
									continue # since it is still the same gene, move on but do not penalize
								else:
									ncontig = refdict[nmatch].scaffold
								if not ncontig==refcontig: # if contigs are different match
									walksteps -= 1 # then skip to next walk step and decrement
									if args.verbose:
										print >> sys.stderr, "# {} matches {} on wrong contig {}, skipping gene".format(ntrans, nmatch, ncontig)
								else:
									if args.verbose:
										print >> sys.stderr, "# Match {} found for {}".format(nmatch, ntrans)
									walksteps = args.span # if a gene is found, reset steps
									endpos = ngeneinfo.end
									accountedgenes.append(ntrans)
									syntenylist.append( (ntrans,nmatch) )
						else: # if walksteps is 0, then break out of for loop
							if args.verbose:
								print >> sys.stderr, "# Limit reached for {}, stopping walk".format(ntrans)
							break # no more searching for genes after walksteps is 0
					### WRITE LONGEST MATCH
					blocklen = len(syntenylist)
					if blocklen >= args.minimum:
						if args.verbose:
							print >> sys.stderr, "# Found block of {0} genes starting from {1}".format(blocklen, startingtrans)
						try:
							if blocklen > max(blocklengths.keys()):
								print >> sys.stderr, "New longest block blk-{} of {} on {}".format(blocknum, blocklen, scaffold)
						except ValueError: # for first check where blocklengths is empty
							print >> sys.stderr, "New longest block blk-{} of {} on {}".format(blocknum, blocklen, scaffold)
						basetotal += transdict[syntenylist[-1][0]].end-transdict[syntenylist[0][0]].start
						blocklengths[blocklen] += 1
						for pair in syntenylist:
							outline = "{}\t{}\tblk-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(scaffold, refcontig, blocknum, pair[0], transdict[pair[0]].start, transdict[pair[0]].end, transdict[pair[0]].strand, pair[1], refdict[pair[1]].start, refdict[pair[1]].end, refdict[pair[1]].strand)
							print >> wayout, outline
						blocknum += 1
					else:
						if args.verbose:
							print >> sys.stderr, "# Block only contained {0} genes, ignoring block".format(blocklen)
				else:
					if args.verbose:
						print >> sys.stderr, "# Too few genes left from {} on scaffold {}, skipping walk".format(startingtrans, scaffold)
				lastmatch = str(refmatch)
	print >> sys.stderr, "# Found {} possible split genes".format(splitgenes), time.asctime()
	print >> sys.stderr, "# Most genes on a query scaffold was {}".format(max(scaffoldgenecounts.keys())), time.asctime() 
	genetotal = sum(x*y for x,y in blocklengths.items())
	print >> sys.stderr, "# Found {} total putative synteny blocks for {} genes".format(sum(blocklengths.values()), genetotal), time.asctime()
	if len(blocklengths)==0: # if no blocks are found
		sys.exit("### NO SYNTENY DETECTED, CHECK GENE ID FORMAT PARAMTERS -Q -D")
	print >> sys.stderr, "# Average block is {:.2f}, longest block was {} genes".format( 1.0*genetotal/sum(blocklengths.values()), max(blocklengths.keys()) ), time.asctime()
	print >> sys.stderr, "# Total block span was {} bases".format(basetotal), time.asctime()

	### MAKE BLOCK HISTOGRAM
	for k,v in sorted(blocklengths.items(),key=lambda x: x[0]):
		print >> sys.stderr, k, v

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
