#!/usr/bin/env python
# microsynteny.py v1.0 2015-10-09
# v1.3 2020-08-18
# v1.4 2022-10-22 bug fix for exon mode
# v1.5 2023-06-14 add presets for genbank GFFs and proteins
# v1.6 2023-07-06 add subject strict option

'''microsynteny.py v1.5 last modified 2023-07-19

microsynteny.py -q query.gtf -d ref_species.gtf -b query_vs_ref_blast.tab -E ../bad_contigs -g -D '_' --blast-query-delimiter '.' > query_vs_ref_microsynteny.tab

    FOR GENERATION OF BLAST DATA TO LINK SETS
blastx -query query.fasta -db ref_prots.fasta -outfmt 6 -evalue 1e-5 > query_vs_ref_blast.tab

    large blast.tab files can be gzipped as .gz

    IN CASE OF NO OUTPUT, CHECK THAT -Q AND -D ARE SET CORRECTLY
    from gtf files, gene_id is extracted
    the blast query names must match the gene_id
    for example, if gene_id is avic.12345 and blast query is avic.12345.1
    set --blast-query-delimiter '.'

    TWO OUTPUT OPTIONS:
    tab delimited text file of 12 columns, consisting of:
query-scaffold  ref-scaffold  block-number
  query-gene  start  end  strand  ref-gene  start  end  strand  score

    OR GFF format file (use --make-gff flag)
    for Parent feature:
    score is length of the block, target is scaffold of subject
    for match_part feature:
    score is bitscore of the blast hit

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

def attributes_to_dict(attributes):
	'''convert GFF attribute string into dictionary of key-value pairs'''
	attrd = {}
	if attributes.find("ID=")>-1 or attributes.find("Parent=")>-1: # indicates GFF3 format
		# if one of the terms does not have = sign, perhaps Note, then ignore
		attrd = dict([(field.strip().split("=",1)) for field in attributes.split(";") if field.count("=")])
	else: # assume GTF format
		try:
			attrd = dict([(field.strip().split(" ",1)) for field in attributes.split(";") if field])
		except ValueError: # catch for Ensembl genomes, which use = but not ID
			attrlist = [field for field in attributes.split(";") if field]
			for attr in attrlist:
				try:
					if attr.count("=")>0:
						attrd.update(dict([attr.strip().split("=")]))
					elif attr.count(" ")>0:
						attrd.update(dict([attr.strip().split(" ")]))
					else: # apparently the field is not delimited
						attrd["NULL"] = attr
				except ValueError: # for in line comments like some Broad Institute gtfs
					# '# At least one base has a quality score < 10'
					sys.stderr.write("WARNING: UNKNOWN ATTRIBUTE: {}\n".format(attr) )
	return attrd

def make_exclude_dict(excludefile):
	'''read file of list of contigs, and return a dict where keys are contig names to exclude'''
	sys.stderr.write("# Reading exclusion list {}  {}\n".format(excludefile, time.asctime() ) )
	exclusion_dict = {}
	for term in open(excludefile,'r'):
		term = term.strip()
		if term[0] == ">":
			term = term[1:]
		exclusion_dict[term] = True
	sys.stderr.write("# Found {} contigs to exclude  {}\n".format( len(exclusion_dict), time.asctime() ) )
	return exclusion_dict

def parse_gtf(gtffile, exons_to_genes, cds_to_genes, excludedict, delimiter, is_genbank, isref=False):
	'''from a gtf, return a dict of dicts where keys are scaffold names, then gene names, and values are gene info as a tuple of start end and strand direction'''
	if gtffile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing {} as gzipped  {}\n".format(gtffile, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing {}  {}\n".format(gtffile, time.asctime() ) )

	if isref: # meaning is db/subject, thus get normal dictionary
		genesbyscaffold = {}
	else:
		genesbyscaffold = defaultdict(dict) # scaffolds as key, then gene name, then genemapping tuple

	feature_counts = defaultdict(int) # count freq of each feature type

	# used if exons_to_genes or cds_to_genes, otherwise should stay empty
	nametoscaffold = {} # in order to get transcript boundaries, store names to scaffolds
	nametostrand = {} # store strand by gene ID
	exonboundaries = defaultdict(list) # make list of tuples of exons by transcript, to determine genes

	# begin parsing
	for line in opentype(gtffile,'rt'):
		line = line.strip()
		if line and not line[0]=="#": # ignore empty lines and comments
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			if excludedict and excludedict.get(scaffold, False):
				continue # skip anything that hits to excludable scaffolds
			feature = lsplits[2]
			feature_start = int(lsplits[3])
			feature_end = int(lsplits[4])

			# count frequency of all features
			feature_counts[feature] += 1

			attributes = lsplits[8]
			attrd = attributes_to_dict(attributes)
			# for top-level features, get bounds directly
			if feature=="transcript" or feature=="mRNA" or feature=="gene":
				if is_genbank: # just use CDS features, that match the protein IDs
					continue
				raw_geneid = attrd.get("ID",None)
				if raw_geneid is None: # try other format
					raw_geneid = attrd.get("gene_id", None)
				if raw_geneid is None:
					print( "ERROR: cannot extract ID= for {}\n{}".format(attributes, line) , file=sys.stderr )
				# if a delimiter is given for either query or db, then split
				if delimiter:
					geneid = geneid.rsplit(delimiter,1)[0]
				else:
					geneid = raw_geneid

				# generate tuple differently for query and db
				if isref: # for db
					refbounds = refgene(scaffold=lsplits[0], start=feature_start, end=feature_end, strand=lsplits[6] )
					genesbyscaffold[geneid] = refbounds
				else: # for query
					boundstrand = querygene(start=feature_start, end=feature_end, strand=lsplits[6] )
					genesbyscaffold[scaffold][geneid] = boundstrand
			# if using exons only, then start collecting exons
			elif (exons_to_genes and feature=="exon"):
				if is_genbank: # just use CDS features, that match the protein IDs
					continue
				raw_geneid = attrd.get("gene_id",None)
				if raw_geneid is None: # try other format
					raw_geneid = attrd.get("name", None)
				# if a delimiter is given for either query or db, then split
				if delimiter:
					geneid = geneid.rsplit(delimiter,1)[0]
				else:
					geneid = raw_geneid
				nametoscaffold[geneid] = scaffold
				nametostrand[geneid] = lsplits[6]
				feature_bounds = ( feature_start , feature_end )
				exonboundaries[geneid].append(feature_bounds) # for calculating gene boundaries
			# if using only CDS
			elif (cds_to_genes and feature=="CDS"):
				if is_genbank: # just use CDS features, that match the protein IDs
					index_id = attrd.get("protein_id",None)
					if index_id is None: # something went wrong, maybe mismatch format
						index_id = attrd.get("Parent",None)
				else:
					index_id = attrd.get("Parent",None)
				nametoscaffold[index_id] = scaffold
				nametostrand[index_id] = lsplits[6]
				feature_bounds = ( feature_start , feature_end )
				exonboundaries[index_id].append(feature_bounds) # for calculating gene boundaries

	# print overall counts of features
	for k in sorted(feature_counts.keys()):
		sys.stderr.write("{}\t{}\n".format(k, feature_counts[k]) )

	# flag for no exons
	if exons_to_genes and len(exonboundaries)==0:
		if feature_counts.get("CDS",0) > 0:
			sys.stderr.write("WARNING: NO EXONS FOUND, {} CDS features detected, try to rerun with -c\n".format( feature_counts.get("CDS") ) )

	# after parsing file, calculate stats and return
	if len(genesbyscaffold) > 0: # even if no-genes was set, this should be more than 0 if genes were in one gtf
		if isref:
			genecount = len(genesbyscaffold)
			sys.stderr.write("# Found {} top-level features (genes/transcripts)  {}\n".format( genecount, time.asctime() ) )
		else:
			# count up number of genes on each scaffold, by length of each value
			# then convert to list, then sum again
			genecount = sum( list( map( len,genesbyscaffold.values() ) ) )
			sys.stderr.write("# Found {} top-level features (genes/transcripts)  {}\n".format( genecount, time.asctime() ) )
		return genesbyscaffold
	# for exon or CDS mode
	else: # generate gene boundaries by scaffold
		subpart_total = sum(len(x) for x in exonboundaries.values())
		sys.stderr.write("# Estimated {} genes from {} subparts  {}\n".format(len(exonboundaries), subpart_total, time.asctime() ) )
		if isref: # make different tuple for reference genes, since they are indexed by name, not scaffold
			for gene,exons in exonboundaries.items():
				refbounds = refgene(scaffold=nametoscaffold[gene], start=min(x[0] for x in exons), end=max(x[1] for x in exons), strand=nametostrand[gene] )
				genesbyscaffold[gene] = refbounds
			genecount = len(genesbyscaffold) # uses len here
		else:
			for gene,exons in exonboundaries.items():
				boundstrand = querygene(start=min(x[0] for x in exons), end=max(x[1] for x in exons), strand=nametostrand[gene] )
				genesbyscaffold[nametoscaffold[gene]][gene] = boundstrand
			genecount = sum(len(x) for x in genesbyscaffold.values()) # must get length for each sub-dict of each scaffold
		sys.stderr.write("# Inferred bounds for {} genes  {}\n".format( genecount, time.asctime() ) )
		return genesbyscaffold

def parse_tabular_blast(blasttabfile, evaluecutoff, querydelimiter, refdelimiter, switchquery, maxhits):
	'''read tabular blast file, return a dict where key is query ID and value is subject ID'''
	if blasttabfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing tabular blast output {} as gzipped  {}\n".format(blasttabfile, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing tabular blast output {}  {}\n".format(blasttabfile, time.asctime() ) )
	query_to_sub_dict = defaultdict(dict)
	sub_counts_dict = defaultdict(int) # key is subject ID, value is count
	query_hits = defaultdict(int) # counter of hits
	evalueRemovals = 0
	for line in opentype(blasttabfile, 'rt'):
		line = line.strip()
		lsplits = line.split("\t")
		# qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
		if switchquery:
			queryseq = lsplits[1].rsplit(refdelimiter,1)[0]
			subjectid = lsplits[0].rsplit(querydelimiter,1)[0]
		else:
			queryseq = lsplits[0].rsplit(querydelimiter,1)[0]
			subjectid = lsplits[1].rsplit(refdelimiter,1)[0]

		# filter by evalue
		if float(lsplits[10]) > evaluecutoff:
			evalueRemovals += 1
			continue
		# filter by number of hits
		if query_hits.get(queryseq,0) >= maxhits: # too many hits already, skip
			continue
		# otherwise add the entry
		bitscore = float(lsplits[11])
		query_to_sub_dict[queryseq][subjectid] = bitscore
		sub_counts_dict[subjectid] += 1
		query_hits[queryseq] += 1
	sys.stderr.write("# Found blast hits for {} query sequences and {} subjects  {}\n".format( len(query_to_sub_dict), len(sub_counts_dict), time.asctime() ) )
	sys.stderr.write("# Removed {} hits by evalue, kept {} hits\n".format( evalueRemovals, sum(query_hits.values()) ) )
	sys.stderr.write("# Names parsed as {} from {}, and {} from {}\n".format( queryseq,lsplits[0], subjectid,lsplits[1] ) )

	# filter by number of hits
	total_kept = 0
	large_group_removals_qu = {} # to prevent multiple counting, store keys
	large_group_removals_sb = {} # or possibly to later check what was removed
	filtered_hit_dict = defaultdict( dict )
	for queryseq, subdict in query_to_sub_dict.items():
		num_hits = len(subdict)
		if num_hits >= maxhits: # remove proteins with many hits, as large protein families likely lead to spurious synteny
			large_group_removals_qu[queryseq] = True
		hit_counter = 0 # reset for each query, to take no more than maxhits
		for subseq, bits in sorted(subdict.items(), key=lambda x: x[1], reverse=True):
			if sub_counts_dict.get(subseq, 0) >= maxhits:
				large_group_removals_sb[subseq] = True
			if hit_counter >= maxhits:
				continue
			if large_group_removals_qu.get(queryseq,False) or large_group_removals_sb.get(subseq,False):
				continue
			filtered_hit_dict[queryseq][subseq] = bits
			hit_counter += 1 # should never get above maxhits
		total_kept += hit_counter
	sys.stderr.write("# Removed {} queries and {} subjects with {} or more hits\n".format( len(large_group_removals_qu), len(large_group_removals_sb), maxhits ) )
	return filtered_hit_dict

def randomize_genes(refdict):
	'''take the gtf dict and randomize the gene names for all genes, return a similar dict of dicts'''
	genepositions = {} # store gene positions as tuples
	randomgenelist = []
	sys.stderr.write("# Randomizing query gene positions  {}\n".format( time.asctime() ) )
	for scaffold, genedict in refdict.items(): # iterate first to get list of all genes
		for genename in genedict.keys():
			randomgenelist.append(genename)
			genepositions[genename] = genedict[genename]
	# randomize the list
	random.shuffle(randomgenelist)
	# reiterate in same order, but store random gene names at the same position
	genecounter = 0
	randomgenesbyscaf = defaultdict(dict) # scaffolds as key, then gene name, then genemapping tuple
	for scaffold, genedict in refdict.items(): # iterate again to reassign genes to each scaffold
		for genename, bounds in genedict.items():
			randomgenesbyscaf[scaffold][randomgenelist[genecounter]] = genepositions[genename]
			genecounter += 1
	sys.stderr.write("# Randomized {} genes  {}\n".format(genecounter, time.asctime() ) )
	return randomgenesbyscaf

def synteny_walk(querydict, blastdict, refdict, min_block, max_span, max_distance, is_subject_strict, is_verbose, wayout, make_gff, is_w ):
	'''for each query scaffold, begin with the first gene and follow as far as possible to identify colinear blocks, then print to stdout'''
	syntenylist = [] # list to keep track of matches, as tuples of (querygene_name, subject_name)
	blocknum = 1
	blocklengths = defaultdict(int) # dictionary to keep track of number of blocks of length N
	scaffoldgenecounts = defaultdict(int)
	matched_query_genes = defaultdict(int) # key is gene ID, value is count of matches
	matched_subject_genes = defaultdict(int) # key is gene ID, value is count of matches
	splitgenes = 0 # counter for number of genes where next gene hits same reference, so query might be split
	basetotal = 0 # counter for block length in bases

	lastmatch = "" # string of name of last gene that matched

	querypos = (0,1) # position of first gene on query scaffold
	subjectpos = (0,1) # position of first gene on ref scaffold

	# for each scaffold
	# iterate through the ordered list of transcripts
	sys.stderr.write("# searching for colinear blocks of at least {} genes, with up to {} intervening genes\n".format( min_block, max_span ) )
	for scaffold, transdict in sorted(querydict.items(), key=lambda x: x[0]):
		orderedtranslist = sorted(transdict.items(), key=lambda x: x[1].start) # sort by start position
		genesonscaff = len(orderedtranslist) # keeping track of number of genes by scaffold, for scale
		scaffoldgenecounts[genesonscaff] += 1
		accounted_query_genes = [] # list of genes already in synteny blocks on this contig
		# accounted_subject_genes is declared later
		if is_verbose:
			sys.stderr.write("#1 Scanning query scaffold {0} with {1} genes\n".format(scaffold, genesonscaff) )
		if genesonscaff < min_block: # not enough genes, thus no synteny would be found
			if is_verbose:
				sys.stderr.write("#1 Only {1} genes on scaffold {0}, skipping scaffold\n".format(scaffold, genesonscaff) )
			continue
		# for each transcript
		# if there are blast matches, more than one, perform a synteny walk for each match
		# this will only keep ones that match the same scaffold at other steps
		# this can be multiple hits on the same contig
		# that is, assuming that a tandem duplication of multiple genes can be found
		# such that genes X Y and Z can blast to A B and C, but also downstream to A' B' and C'
		for i, querygene_tuple in enumerate(orderedtranslist):
			# starting from each transcript
			querygene = querygene_tuple[0]
			startingtrans = querygene
			if is_verbose:
				print("## gene {}  prot {}  scaffold {}".format(querygene, startingtrans, scaffold ), file=sys.stderr)
			querypos = (querygene_tuple[1].start, querygene_tuple[1].end)
			# get the blast matches of the query
			blastrefmatch_dict = blastdict.get(startingtrans,None) 
			if blastrefmatch_dict==None: # if no blast match, then skip to next gene
				if is_verbose:
					sys.stderr.write("#2 No blast matches for {}, skipping walk\n".format(startingtrans) )
				continue
			# otherwise start iterating through all blast hits of that gene
			for blast_refgene, bitscore1 in sorted(blastrefmatch_dict.items(), key=lambda x: x[1], reverse=True):
				walksteps = max_span # genes until drop, walk starts new for each transcript
				# get scaffold and position of matched gene
				refscaffold = refdict[blast_refgene].scaffold
				subjectpos = (refdict[blast_refgene].start, refdict[blast_refgene].end)
				if blast_refgene==lastmatch: # unique test to see if blast hit matches previous
					if is_verbose:
						sys.stderr.write("#2 Gene {} matches last gene {}, maybe long split  s{}\n".format(startingtrans, lastmatch, walksteps) )
					splitgenes += 1

				# if gene is already in a synteny block, ignore it
				# due to multiple hits within the same protein, e.g. multidomain proteins
				if startingtrans in accounted_query_genes:
					if is_verbose:
						sys.stderr.write("#2 Gene {} already used in block on {}, skipping  s{}\n".format(startingtrans, scaffold, walksteps) )
					continue
				if startingtrans in matched_query_genes:
					if is_subject_strict: # do not allow double matches to ANY subject
						if is_verbose:
							sys.stderr.write("#2 Gene {} in another block on {}, skipping  s{}\n".format(startingtrans, scaffold, walksteps) )
						continue
				if blast_refgene in matched_subject_genes:
					if is_subject_strict: # do not allow double matches to ANY subject
						if is_verbose:
							sys.stderr.write("#2 Match {} in another block on {}, skipping  s{}\n".format(blast_refgene, refscaffold, walksteps) )
						continue

				# renew synteny list for each query gene
				syntenylist = [ (startingtrans,blast_refgene) ]
				accounted_subject_genes = [blast_refgene]
				prev_match = blast_refgene

				if i < genesonscaff - 1: # this allows for 2 genes left
					######################################
					# begin of gene walk on forward strand
					if is_verbose:
						sys.stderr.write("#3 Starting walk from gene {} on scaffold {} against {}  s{}\n".format(startingtrans, scaffold, blast_refgene, walksteps) )
					for next_gene, next_gene_info in orderedtranslist[i+1:]: # getting next transcript, and next gene
						if walksteps <= 0: # if walksteps is 0, then break out of for loop
							if is_verbose:
								sys.stderr.write("#3 Limit reached at {}, stopping walk for {}\n".format(next_gene, startingtrans) )
							break # no more searching for genes after walksteps is 0
						dist_to_next_query = next_gene_info.start-querypos[1]
						if dist_to_next_query > max_distance: # next gene is too far
							if is_verbose:
								sys.stderr.write("#3 Next query gene {} bp away from {}, stopping walk\n".format(dist_to_next_query, next_gene) )
							break # end gene block
						# update query positions
						querypos = (next_gene_info.start, next_gene_info.end)
						next_match_dict = blastdict.get(next_gene,None)
						# if no blast match, then skip to next walk step, and decrement
						if next_match_dict==None:
							walksteps -= 1
							if is_verbose:
								sys.stderr.write("#3 No blast matches for {}, skipping gene  s{}\n".format(next_gene, walksteps) )
							continue
						# otherwise iterate through matches, finding one within range
						for next_match, bitscoreN in sorted(next_match_dict.items(), key=lambda x: x[1], reverse=True):
							if next_match==prev_match: # if last match is the same as current match, meaning split query
								splitgenes += 1
								accounted_query_genes.append(next_gene)
								if is_verbose:
									sys.stderr.write("#4 {} match {} in the same block on {}, skipping  s{}\n".format( next_gene, next_match, scaffold, walksteps) )
								continue # since it is still the same gene, move on but do not penalize
							else:
								next_ref_scaf = refdict[next_match].scaffold
								next_ref_pos = (refdict[next_match].start, refdict[next_match].end)
							if next_ref_scaf==refscaffold: # scaffolds match
								if next_match in accounted_subject_genes:
									if is_subject_strict: # do not allow double matches to ANY subject
										if is_verbose:
											sys.stderr.write("#4 {} match {} in the same block on {}, skipping  s{}\n".format( next_gene, next_match, scaffold, walksteps) )
										continue
								# determine distance, strand is not considered
								# meaning one value should be negative, other should be positive
								dist_to_next_ref = next_ref_pos[0] - subjectpos[1]
								dist_to_prev_ref = subjectpos[0] - next_ref_pos[1]
								# if either are greater than max distance, then gene is too far
								if dist_to_next_ref > max_distance or dist_to_prev_ref > max_distance:
									if is_verbose:
										sys.stderr.write("#4 {} match to {} is too far, {}bp, ignoring match\n".format(next_gene, next_match, max([dist_to_next_ref,dist_to_prev_ref]) ) )
									continue
								if is_verbose:
									sys.stderr.write("#5 Match {} found for {} on {}  s{}\n".format(next_match, next_gene, scaffold, walksteps ) )
								walksteps = max_span # if a gene is found, reset steps
								subjectpos = next_ref_pos
								prev_match = next_match
								accounted_query_genes.append(next_gene)
								accounted_subject_genes.append(next_match)
								syntenylist.append( (next_gene,next_match) )
								break # move to next gene
							else:
								if is_verbose:
									sys.stderr.write("#4 {} matches {} on wrong contig {}, skipping gene\n".format(next_gene, next_match, next_ref_scaf) )
						else: # ends next_match for loop, meaning no match was found on correct contig
							walksteps -= 1 # then skip to next walk step and decrement

					# when cannot loop further
					### WRITE LONGEST MATCH
					blocklen = len(syntenylist)
					if blocklen >= min_block:
						if is_verbose:
							sys.stderr.write("# Found block blk-{3} of {0} genes starting from {1} on {2}\n".format(blocklen, startingtrans, scaffold, blocknum) )
						try:
							if blocklen > max(blocklengths.keys()):
								sys.stderr.write("New longest block blk-{} of {} on {}\n".format(blocknum, blocklen, scaffold) )
						except ValueError: # for first check where blocklengths is empty
							sys.stderr.write("New longest block blk-{} of {} on {}\n".format(blocknum, blocklen, scaffold) )
						qblockstart = transdict[syntenylist[0][0]].start
						qblockend = transdict[syntenylist[-1][0]].end
						sblockstart = refdict[syntenylist[0][1]].start
						sblockend = refdict[syntenylist[-1][1]].end
						basetotal += qblockend - qblockstart
						blocklengths[blocklen] += 1
						if sblockend > sblockstart:
							strand = "+"
						else:
							strand = "-"
							sblockstart, sblockend = sblockend, sblockstart
						###############################
						# generate GFF for entire block
						if make_gff:
							# could also be "cross_genome_match"
							blockline = "{0}\tmicrosynteny\tmatch\t{1}\t{2}\t{3}\t{4}\t.\tID=blk-{5};Name=blk-{5}_to_{6};Target={6} {7} {8}\n".format( scaffold, qblockstart, qblockend, blocklen, strand, blocknum, refscaffold, sblockstart, sblockend)
							wayout.write(blockline)
							# make GFF line for each match
							for j,pair in enumerate(syntenylist):
								matched_query_genes[pair[0]] += 1
								matched_subject_genes[pair[1]] += 1
								outline = "{0}\tmicrosynteny\tmatch_part\t{1}\t{2}\t{3}\t{4}\t.\tID=blk-{5}.{10}.{11};Parent=blk-{5};Target={6} {7} {8} {9}\n".format( scaffold, transdict[pair[0]].start, transdict[pair[0]].end, blastdict[pair[0]][pair[1]], transdict[pair[0]].strand, blocknum, pair[1], refdict[pair[1]].start, refdict[pair[1]].end, refdict[pair[1]].strand, j+1, pair[0])
								wayout.write(outline)
						############################
						# otherwise use output of v1
						else:
							for pair in syntenylist:
								matched_query_genes[pair[0]] += 1
								matched_subject_genes[pair[1]] += 1
								outline = "{}\t{}\tblk-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(scaffold, refscaffold, blocknum, pair[0], transdict[pair[0]].start, transdict[pair[0]].end, transdict[pair[0]].strand, pair[1], refdict[pair[1]].start, refdict[pair[1]].end, refdict[pair[1]].strand, blastdict[pair[0]][pair[1]])
								wayout.write(outline)
						blocknum += 1
					else:
						if is_verbose:
							sys.stderr.write("# Block only contained {0} genes, ignoring block\n".format(blocklen) )
				else:
					if is_verbose:
						sys.stderr.write("# Too few genes left from {} on scaffold {}, skipping walk\n".format(startingtrans, scaffold) )
				lastmatch = str(blast_refgene)
	final_block_count = sum(list(blocklengths.values()))
	sys.stderr.write("# Found {} possible split genes  {}\n".format(splitgenes, time.asctime() ) )
	sys.stderr.write("# Most genes on a query scaffold was {}\n".format(max(list(scaffoldgenecounts.keys() ) ) ) )
	genetotal = sum(x*y for x,y in blocklengths.items())
	sys.stderr.write("# Found {} total putative synteny blocks for {} genes (may include duplicates)\n".format(final_block_count, genetotal) )
	sys.stderr.write("# Included {} queries and {} target genes\n".format( len(matched_query_genes), len(matched_subject_genes) ) )
	if len(blocklengths)==0: # if no blocks are found
		sys.exit("### NO SYNTENY DETECTED, CHECK GENE ID FORMAT PARAMTERS -Q -D")
	sys.stderr.write("# Average block is {:.2f}, longest block was {} genes\n".format( 1.0*genetotal/final_block_count, max(list(blocklengths.keys())) ) )
	sys.stderr.write("# Total block span was {} bases  {}\n".format(basetotal, time.asctime() ) )

	# w mode output
	if is_w: # summarize on one line
		sys.stderr.write( "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n".format( max_distance, min_block, max_span, final_block_count, sum([len(v) for v in querydict.values()]), len(refdict), len(matched_query_genes), len(matched_subject_genes), 1.0*genetotal/final_block_count, max(list(blocklengths.keys())) ) )
	### MAKE BLOCK HISTOGRAM
	else:
		for k,v in sorted(blocklengths.items(),key=lambda x: x[0]):
			sys.stderr.write("{}\t{}\n".format(k, v) )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast', help="tabular blast output", required=True)
	parser.add_argument('-q','--query-gtf', help="gtf file of query genes", required=True)
	parser.add_argument('-d','--db-gtf', help="gtf file of reference genes", required=True)
	parser.add_argument('-Q','--query-delimiter', help="gene transcript separator for query [.]")
	parser.add_argument('-D','--db-delimiter', help="gene transcript separator for db [.]")
	parser.add_argument('--blast-query-delimiter', help="gene transcript separator for blast query [|]", default='|')
	parser.add_argument('--blast-db-delimiter', help="gene transcript separator for blast ref [|]", default='|')
	parser.add_argument('-c','--cds-only', action="store_true", help="genes are not defined, get gene ID for each CDS feature")
	parser.add_argument('-e','--evalue', type=float, default=1e-4, help="evalue cutoff for post blast filtering [1e-4]")
	parser.add_argument('-E','--exclude', help="file of list of bad contigs, from either genome")
	parser.add_argument('-g','--no-genes', action="store_true", help="genes are not defined, get gene ID for each exon")
	parser.add_argument('-m','--minimum', type=int, default=3, help="minimum syntenic genes to keep block, must be >=2 [3]")
	parser.add_argument('-G','--group-size-maximum', metavar="N", type=int, default=100, help="remove queries with more than N hits, e.g. transposons [100]")
	parser.add_argument('-s','--span', type=int, default=5, help="max number of skippable genes [5]")
	parser.add_argument('-z','--distance', type=int, default=30000, help="max distance on query scaffold before next gene [30000]")
	parser.add_argument('--make-gff', help="make GFF output, instead of tabular blocks", action="store_true")
	parser.add_argument('--genbank-gff', help="use presets when proteins and GFF files are from GenBank", action="store_true")
	parser.add_argument('-R','--randomize', help="randomize positions of query GTF", action="store_true")
	parser.add_argument('-S','--switch-query', help="switch query and subject", action="store_true")
	parser.add_argument('-T','--strict', help="require strict synteny", action="store_true")
	parser.add_argument('-v','--verbose', help="verbose output", action="store_true")
	parser.add_argument('-w','--w', action="store_true", help="give summarized output")
	args = parser.parse_args(argv)

	sys.stderr.write("# Running command:\n{}\n".format( ' '.join(sys.argv) ) )

	if args.minimum < 2:
		sys.stderr.write("WARNING: MINIMUM COLINEARITY -m MUST BE GREATER THAN 1, {} GIVEN\n".format(args.minimum) )
		sys.stderr.write("SETTING MINIMUM COLINEARITY TO 2\n")
		args.minimum = 2

	exclusiondict = make_exclude_dict(args.exclude) if args.exclude else {}

	# note user settings
	if args.cds_only:
		sys.stderr.write("# -c ENABLED, will determine genes from CDS features\n")
	if args.genbank_gff:
		args.cds_only = True
	### SETUP DICTIONARIES ###
	if args.switch_query:
		querydict = parse_gtf(args.db_gtf, args.no_genes, args.cds_only, exclusiondict, args.db_delimiter, args.genbank_gff, isref=False)
		refdict = parse_gtf(args.query_gtf, args.no_genes, args.cds_only, exclusiondict, args.query_delimiter, args.genbank_gff, isref=True)
		blastdict = parse_tabular_blast(args.blast, args.evalue, args.blast_query_delimiter, args.blast_db_delimiter, args.switch_query, args.group_size_maximum)
	else:
		querydict = parse_gtf(args.query_gtf, args.no_genes, args.cds_only, exclusiondict, args.query_delimiter, args.genbank_gff, isref=False )
		refdict = parse_gtf(args.db_gtf, args.no_genes, args.cds_only, exclusiondict, args.db_delimiter, args.genbank_gff, isref=True)
		blastdict = parse_tabular_blast(args.blast, args.evalue, args.blast_query_delimiter, args.blast_db_delimiter, args.switch_query, args.group_size_maximum )

	### IF DOING RANDOMIZATION ###
	if args.randomize:
		querydict = randomize_genes(querydict)

	if args.make_gff:
		sys.stderr.write("# make GFF output: {}\n".format( args.make_gff ) )
	### START SYNTENY WALKING ###
	synteny_walk(querydict, blastdict, refdict, args.minimum, args.span, args.distance, args.strict, args.verbose, wayout, args.make_gff, args.w )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
