#!/usr/bin/env python
#
# scaffold_synteny.py created 2019-03-27

'''
scaffold_synteny.py  v1.1 last modified 2019-10-08
    makes a table of gene matches between two genomes, to detect synteny
    these can be converted into a dotplot of gene matches

    options are:
    -b : tabular blast output
    -f and -F : fasta files of scaffolds for query and subject genomes
    -q and -d : GFF files of genes for query and subject genomes
    -l and -L : total length to keep for each genome, in MB
        i.e. to only take long scaffolds up to length -l

scaffold_synteny.py -b monbr1_vs_srosetta_blastp.tab -q Monbr1_augustus_v1_no_comment.gff -d Srosetta_mrna_only_ID_renamed.gff -f Monbr1_scaffolds.fasta -F Salpingoeca_rosetta.dna.toplevel.fa.gz -l 40 -L 50 > monbr1_vs_srosetta_scaffold2d_points.tab

    output consists of 8 columns, which differ by data, but are all
      in the same output file
    first type is scaffold information, for both genomes, as:
  genome (either 1 or 2)  scaffold-name  number  length
    fraction of total assembly  cumulative sum  cumulative fraction  blank

s1  contig_001_length_2062630  1  2062630  0.023656  2062630  0.023656  -

    second type is positions of matches, relative to total length:
  symbol indicating gene  gene  contig of gene  match to other genome
    scaffold of match  position in genome 1  position in genome 2  bitscore

g  braker1_g09939  contig_176_length_141897  g5612.t1  scaffold_5  72253822  38899018  104.0

   generate tabular blast data with blastp:
blastp -query Monbr1_augustus_v1.prot_no_rename.fasta -db Salpingoeca_rosetta.pep.all.fa -outfmt 6 -num_threads 6 -evalue 1e-3 -max_target_seqs 100 > monbr1_vs_srosetta_blastp.tab

   generate plot using synteny_2d_plot.R:
Rscript synteny_2d_plot.R monbr1_vs_srosetta_scaffold2d_points.tab Monosiga-brevicollis Salpingoeca-rosetta
'''

import sys
import argparse
import time
import os
import re
import gzip
import random
from collections import defaultdict
from Bio import SeqIO

def make_seq_length_dict(contigsfile, maxlength, exclusiondict, wayout, isref=False):
	'''read fasta file, and return dict where key is scaffold name and value is length'''
	lengthdict = {}
	if contigsfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing genomic contigs {} as gzipped  ".format(contigsfile) + time.asctime() + os.linesep)
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing genomic contigs {}  ".format(contigsfile) + time.asctime() + os.linesep)
	for seqrec in SeqIO.parse(opentype(contigsfile,'rt'), "fasta"):
		lengthdict[seqrec.id] = len(seqrec.seq)
	sys.stderr.write("# Found {} contigs  ".format(len(lengthdict)) + time.asctime() + os.linesep)

	# make scaffold key as s1 for query and s2 for reference db
	scafkey = "s2" if isref else "s1"
	totalgenomesize = sum( list(lengthdict.values() ) )

	sorteddict = {} # key is scaffold, value is offset relative to previous scaffolds
	lengthsum = 0 # cumulative sum of lengths of large scaffolds
	breaklines = [0] # contains scaffold lengths, starts with position 0
	maxlength_MB = 1000000*maxlength
	sys.stderr.write("# Sorting contigs by length, keeping up to {}Mbp\n".format(maxlength) )
	# keep the first N scaffolds, where total length is approximately maxlength
	scaffoldcounter = 0
	for k,v in sorted(lengthdict.items(), key=lambda x: x[1], reverse=True):
		if k in exclusiondict:
			continue
		scaffoldcounter += 1
		sorteddict[k] = lengthsum
		lengthsum += v
		wayout.write("{}\t{}\t{}\t{}\t{:.6f}\t{}\t{:.6f}\t-\n".format(scafkey, k, scaffoldcounter, v, v*1.0/totalgenomesize, lengthsum, lengthsum*1.0/totalgenomesize) )
		if lengthsum >= maxlength_MB:
			break
	#	breaklines.append(lengthsum)
	# offset len(sorteddict) count by 1, as length 0 is not a contig
	sys.stderr.write("# Kept {} contigs, for {} bases, last contig was {}bp long  ".format( len(sorteddict)-1, lengthsum, v) + time.asctime() + os.linesep)
	#sys.stderr.write(",".join( map(str,breaklines) )	
	return sorteddict

def parse_gtf(gtffile, excludedict, delimiter, isref=False):
	'''from a gtf, return a dict of dicts where keys are scaffold names, then gene names, and values are gene info as a tuple of start end and strand direction'''
	if isref:
		genesbyscaffold = {} # key is reference gene ID, value is list of scaffold and gene midpoint position
	else:
		genesbyscaffold = defaultdict(dict) # scaffolds as key, then gene name, then gene position integer

	if gtffile.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing loci from {} as gzipped  ".format(gtffile) + time.asctime() + os.linesep)
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing loci from {}  ".format(gtffile) + time.asctime() + os.linesep)
	for line in opentype(gtffile,'rt'):
		line = line.strip()
		if line and not line[0]=="#": # ignore empty lines and comments
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			if excludedict and excludedict.get(scaffold, False):
				continue # skip anything that hits to excludable scaffolds
			feature = lsplits[2]
			attributes = lsplits[8]
			if feature=="gene" or feature=="transcript" or feature=="mRNA":
				geneid = re.search('ID=([\w.|-]+)', attributes).group(1)
				# if a delimiter is given for either query or db, then split
				if delimiter:
					geneid = geneid.rsplit(delimiter,1)[0]

				# generate midpoint of each gene as average of start and end positions
				genemidpoint = (int(lsplits[3]) + int(lsplits[4])) / 2
				if isref:
					genesbyscaffold[geneid] = [scaffold,genemidpoint]
				else:
					genesbyscaffold[scaffold][geneid] = genemidpoint

	if len(genesbyscaffold) > 0:
		sys.stderr.write("# Found {} genes  ".format( sum( list( map( len,genesbyscaffold.values()) ) ) ) + time.asctime() + os.linesep)
		return genesbyscaffold
	else:
		sys.stderr.write("# WARNING: NO GENES FOUND\n")

def parse_tabular_blast(blasttabfile, evaluecutoff, querydelimiter, refdelimiter, maxhits=1):
	'''read tabular blast file, return a dict where key is query ID and value is dict of subject ID and bitscore'''
	if blasttabfile.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing tabular blast output {} as gzipped  ".format(blasttabfile) + time.asctime() + os.linesep)
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing tabular blast output {}  ".format(blasttabfile) + time.asctime() + os.linesep)
	query_to_sub_dict = defaultdict( lambda: defaultdict(int) )
	query_hits = defaultdict(int) # counter of hits
	evalueRemovals = 0
	for line in opentype(blasttabfile, 'rt'):
		line = line.strip()
		lsplits = line.split("\t")
		# qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
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
		query_to_sub_dict[queryseq][subjectid] += bitscore
		query_hits[queryseq] += 1
	sys.stderr.write("# Found blast hits for {} query sequences  ".format( len(query_to_sub_dict) ) + time.asctime() + os.linesep)
	sys.stderr.write("# Removed {} hits by evalue, kept {} hits".format( evalueRemovals, sum( list(query_hits.values())) ) + os.linesep)
	sys.stderr.write("# Names parsed as {} from {}, and {} from {}".format( queryseq,lsplits[0], subjectid,lsplits[1] ) + os.linesep)
	return query_to_sub_dict

def generate_synteny_points(queryScafOffset, dbScafOffset, queryPos, dbPos, blastdict, wayout):
	'''combine all datasets and for each gene on the query scaffolds, print tab delimited data'''
	printcount = 0
	scaffoldtotals = defaultdict(int) # counts of total genes for each scaffold
	sys.stderr.write("# Determining match positions  " + time.asctime() + os.linesep)
	for scaffold, genedict in queryPos.items():
		scaffoldcounts = defaultdict(int) # counts of hits to each reference scaffold
		for gene, localposition in genedict.items():
			scaffoldtotals[scaffold] += 1
			queryoffset = queryScafOffset.get(scaffold,None)
			if queryoffset is None:
				continue
			overallposition = localposition + queryoffset
			blasthits = blastdict.get(gene, None)
			if blasthits is None:
				continue
			for matchgene, bitscore in sorted(blasthits.items(), key=lambda x: x[1], reverse=True):
				matchscaf, matchposition = dbPos.get(matchgene, [None, None])
				if matchscaf is None:
					continue
				matchoffset = dbScafOffset.get(matchscaf, None)
				if matchoffset is None:
					continue
				overallmatchpos = matchposition + matchoffset
				scaffoldcounts[matchscaf] += 1
				printcount += 1
				wayout.write("g\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, scaffold, matchgene, matchscaf, overallposition, overallmatchpos, bitscore) )
	if printcount:
		sys.stderr.write("# Wrote match positions for {} genes\n".format( printcount ) )
	else:
		sys.stderr.write("# WARNING: NO MATCHES FOUND, CHECK -Q AND -D\n")

def randomize_genes_globally(refdict):
	'''take the query gtf dict and randomize the gene names for all genes on all scaffolds, return a similar dict of dicts'''
	genepositions = {} # store gene positions
	randomgenelist = []
	sys.stderr.write("# Globally randomizing query gene positions  " + time.asctime() + os.linesep)
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
	sys.stderr.write("# Randomized {} genes  ".format(genecounter) + time.asctime() + os.linesep)
	return randomgenesbyscaf

def randomize_genes_locally(refdict):
	'''take the query gtf dict and randomize the gene names for genes within each scaffold, return a similar dict of dicts'''
	genepositions = {} # store gene positions
	scaffoldcount = 0
	total_genes = 0
	randomgenesbyscaf = defaultdict(dict) # scaffolds as key, then gene name, then genemapping tuple
	sys.stderr.write("# Randomizing query gene positions by scaffold  " + time.asctime() + os.linesep)
	for scaffold, genedict in refdict.items(): # iterate first to get list of all genes
		scaffoldcount += 1
		randomgenelist = [] # generate new list to randomize for each scaffold
		for genename in genedict.keys():
			randomgenelist.append(genename)
			genepositions[genename] = genedict[genename]
		# randomize the list
		random.shuffle(randomgenelist)
		# reiterate in same order, but store random gene names at the same position
		genecounter = 0
		for genename, bounds in genedict.items(): # iterate again to reassign genes to each scaffold
			randomgenesbyscaf[scaffold][randomgenelist[genecounter]] = genepositions[genename]
			genecounter += 1
			total_genes += 1
	sys.stderr.write("# Randomized {} genes on {} scaffolds  ".format(total_genes,scaffoldcount) + time.asctime() + os.linesep)
	return randomgenesbyscaf

def randomize_db_locally(refdict):
	'''take the ref gtf dict and randomize the gene positions for genes within each scaffold, return a similar dict of lists'''
	scaffoldcount = 0
	total_genes = 0
	randomgenesbyscaf = {} # key is reference gene ID, value is list of scaffold and gene midpoint position
	sys.stderr.write("# Randomizing reference gene positions by scaffold  " + time.asctime() + os.linesep)
	for scaffold, genedict in refdict.items(): # iterate first to get list of all genes
		scaffoldcount += 1
		randomposlist = [] # generate new list of the gene positions to randomize for each scaffold
		for genename in genedict.keys():
			randomposlist.append( genedict[genename] )
		# randomize the list
		random.shuffle(randomposlist)
		# reiterate in same order, but store random positions for each gene
		genecounter = 0
		for genename, bounds in genedict.items(): # iterate again to reassign genes to each scaffold
			randomgenesbyscaf[genename] = [scaffold, randomposlist[genecounter]]
			genecounter += 1
			total_genes += 1
	sys.stderr.write("# Randomized {} genes on {} scaffolds  ".format(total_genes,scaffoldcount) + time.asctime() + os.linesep)
	return randomgenesbyscaf

def make_exclude_dict(excludefile):
	'''read file of list of contigs, and return a dict where keys are contig names to exclude'''
	sys.stderr.write("# Reading exclusion list {}  ".format(excludefile) + time.asctime() + os.linesep)
	exclusion_dict = {}
	for term in open(excludefile,'r'):
		term = term.strip()
		if term[0] == ">":
			term = term[1:]
		exclusion_dict[term] = True
	sys.stderr.write("# Found {} contigs to exclude  ".format(len(exclusion_dict) ) + time.asctime() + os.linesep)
	return exclusion_dict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast', help="tabular blast output", required=True)
	parser.add_argument('-f','--query-fasta', help="fasta file of query scaffolds", required=True)
	parser.add_argument('-F','--db-fasta', help="fasta file of reference scaffolds", required=True)
	parser.add_argument('-q','--query-gff', help="GFF file of query genes", required=True)
	parser.add_argument('-d','--db-gff', help="GFF file of reference genes", required=True)
	parser.add_argument('-Q','--query-delimiter', help="gene transcript separator for query [.]")
	parser.add_argument('-D','--db-delimiter', help="gene transcript separator for db [.]")
	parser.add_argument('--blast-query-delimiter', help="gene transcript separator for blast query [|]", default='|')
	parser.add_argument('--blast-db-delimiter', help="gene transcript separator for blast ref [|]", default='|')
	parser.add_argument('-E','--exclude', help="file of list of bad contigs, from either genome")
	parser.add_argument('-e','--evalue', type=float, default=1e-4, help="evalue cutoff for post blast filtering [1e-4]")
	parser.add_argument('-c','--coverage', default=0.8, type=int, help="minimum alignment coverage [0.8]")
	parser.add_argument('-l','--query-genome-len', type=int, default=100, help="length of query scaffolds, in Mbp [100]")
	parser.add_argument('-L','--db-genome-len', type=int, default=100, help="length of reference scaffolds, in Mbp [100]")
	parser.add_argument('-R','--global-randomize', help="globally randomize gene positions of query GFF, cannot use with -S", action="store_true")
	parser.add_argument('-S','--scaffold-randomize', help="randomize gene positions of query GFF within each scaffold, cannot use with -R", action="store_true")
	parser.add_argument('--double-randomize', help="randomize gene positions of db, use with -S", action="store_true")
	args = parser.parse_args(argv)

	exclusiondict = make_exclude_dict(args.exclude) if args.exclude else {}

	# read both sets of scaffolds
	query_scaf_lengths = make_seq_length_dict(args.query_fasta, args.query_genome_len, exclusiondict, wayout, False)
	db_scaf_lengths = make_seq_length_dict(args.db_fasta, args.db_genome_len, exclusiondict, wayout, True)

	# read query as normal
	query_gene_pos = parse_gtf(args.query_gff, exclusiondict, args.query_delimiter, False)
	### IF DOING RANDOMIZATION ###
	if args.global_randomize:
		query_gene_pos = randomize_genes_globally(query_gene_pos)
	elif args.scaffold_randomize:
		query_gene_pos = randomize_genes_locally(query_gene_pos)

	### IF DOING DOUBLE RANDOMIZE ###
	if args.double_randomize: # read db first as query format, randomize and generate the dict in the ref format
		db_gene_pos = parse_gtf(args.db_gff, exclusiondict, args.db_delimiter, False)
		db_gene_pos = randomize_db_locally(db_gene_pos)
	# if NOT RANDOMIZING REFERENCE #
	else: # otherwise read as normal into the ref format
		db_gene_pos = parse_gtf(args.db_gff, exclusiondict, args.db_delimiter, True)

	# read blast hits
	blastdict = parse_tabular_blast(args.blast, args.evalue, args.blast_query_delimiter, args.blast_db_delimiter)

	# write output
	generate_synteny_points( query_scaf_lengths, db_scaf_lengths, query_gene_pos, db_gene_pos, blastdict, wayout)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
