#!/usr/bin/env python

# blast2genomegff.py v1.0 2016-05-16
#
# and using information from:
# http://www.sequenceontology.org/gff3.shtml

# for SOFA terms:
# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo

'''blast2genomegff.py  last modified 2016-09-29
    convert blast output to gff format for genome annotation
    blastx of a transcriptome (genome guided or de novo) against a protein DB:

blast2genomegff.py -b blastx_out6.tab -d protein_db.fasta -g transcripts.gtf > output.gff

    generate the tabular blastx output (-outfmt 6) by:
blastx -query transcripts.fasta -db protein_db.fasta -outfmt 6 -max_target_seqs 5 > blastx_out6.tab

    change the second and third fields in the gff output with -p and -t

blast2genomegff.py -b blastn_output.tab -p BLASTN -t EST_match > output.gff

    -t should be a standard sequence ontology term, including:
      nucleotide_match   EST_match   protein_match

    BLAST output should be made from blast programs using -outfmt 6

    evalue cutoff between 1e-3 and 1e-40 is appropriate to filter bad hits
    though this depends on the bitscore and so the relatedness of the species

    for AUGUSTUS GFF format genes, use -x and -T
'''

import sys
import argparse
import time
import re
from collections import defaultdict
from Bio import SeqIO

def make_seq_length_dict(sequencefile):
	print >> sys.stderr, "# Parsing target sequences from {}".format(sequencefile), time.asctime()
	lengthdict = {}
	for seqrec in SeqIO.parse(sequencefile,'fasta'):
		lengthdict[seqrec.id] = len(seqrec.seq)
	print >> sys.stderr, "# Found {} sequences".format(len(lengthdict)), time.asctime()
	return lengthdict

def gtf_to_intervals(gtffile, keepcds, transdecoder, nogenemode=False):
	'''convert protein or gene intervals from gff to dictionary where mrna IDs are keys and lists of intervals are values'''
	# this should probably be a class
	geneintervals = defaultdict(list)
	genestrand = {}
	genescaffold = {}

	commentlines = 0
	linecounter = 0
	transcounter = 0
	exoncounter = 0
	print >> sys.stderr, "# Parsing gff from {}".format(gtffile), time.asctime()
	for line in open(gtffile).readlines():
		line = line.strip()
		if line: # ignore empty lines
			if line[0]=="#": # count comment lines, just in case
				commentlines += 1
			else:
				linecounter += 1
				lsplits = line.split("\t")
				scaffold = lsplits[0]
				feature = lsplits[2]
				attributes = lsplits[8]
			#	if aqumode:
			#		attributes = attributes.replace('CDS:','') # remove CDS: from Aqu2 gene models
				if attributes.find("ID")==0: # indicates gff3 format
					geneid = re.search('ID=([\w.|-]+)', attributes).group(1)
				elif attributes.find("Parent")==0: # gff3 format but no ID
					geneid = re.search('Parent=([\w.|-]+)', attributes).group(1)
				elif attributes.find("gene_id")==0: # indicates gtf format
					geneid = re.search('transcript_id "([\w.|-]+)";', attributes).group(1)
			#	elif jgimode and attributes.find("name")==0: # for JGI like: name "fgeneshTA2_pg.C_scaffold_1000001";
			#		geneid = re.search('name "([\w.|-]+)";', attributes).group(1)
				if transdecoder: # meaning CDS IDs will start with cds.gene.123|m.1
					geneid = geneid.replace("cds.","") # simply remove the cds.
					geneid = geneid.replace(".cds","") # also works for AUGUSTUS
				if feature=="transcript" or feature=="mRNA": # or (aqumode and feature=="gene"):
					transcounter += 1
					strand = lsplits[6]
					genestrand[geneid] = strand
					genescaffold[geneid] = scaffold
				elif feature=="exon" or (keepcds and feature=="CDS"):
					exoncounter += 1
					boundaries = ( int(lsplits[3]), int(lsplits[4]) )
					if nogenemode: # gtf contains only exon and CDS, so get gene info from each CDS
						geneid = re.search('transcript_id "([\w.|-]+)";', attributes).group(1)
					geneintervals[geneid].append(boundaries)
	print >> sys.stderr, "# Counted {} lines and {} comments".format(linecounter, commentlines), time.asctime()
	print >> sys.stderr, "# Counted {} exons for {} transcripts".format(exoncounter, transcounter), time.asctime()
	return geneintervals, genestrand, genescaffold

def parse_tabular_blast(blastfile, lengthcutoff, evaluecutoff, bitscutoff, programname, outputtype, donamechop, swissprot, seqlengthdict, geneintervals, genestrand, genescaffold, debugmode=False):
	'''parse blast hits from tabular blast and write to stdout as genome gff'''
	querynamedict = {} # counter of unique queries
	shortRemovals = 0
	evalueRemovals = 0
	bitsRemovals = 0
	intervalproblems = 0 # counter if no intervals are found for some sequence
	intervalcounts = 0
	backframecounts = 0
	# set up parameters by blast program
	blastprogram = programname.lower()
	if blastprogram=="blastn" or blastprogram=="blastx" or blastprogram=="tblastx":
		multiplier = 1
	else: # meaning blastp or tblastn
		multiplier = 3

	hitDictCounter = defaultdict(int)
	linecounter = 0
	print >> sys.stderr, "# Starting BLAST parsing on %s" % (blastfile), time.asctime()
	for line in open(blastfile, 'r').readlines():
		line = line.strip()
		if not line or line[0]=="#": # skip comment lines
			continue # also catch for empty line, which would cause IndexError
		linecounter += 1
		lsplits = line.split()
		#qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
		sseqid = lsplits[1]

		# do all removals before any real counting
		evalue = float(lsplits[10])
		bitscore = float(lsplits[11])
		alignlength = float(lsplits[3])
		hitstart = int(lsplits[6])
		hitend = int(lsplits[7])
		hitlength = abs(hitend - hitstart) + 1 # bases 1 to 6 should have length 6
		fractioncov = alignlength/seqlengthdict[sseqid]
		bitslength = bitscore/alignlength
		if fractioncov < lengthcutoff: # skip domains that are too short
			shortRemovals += 1
			continue
		if bitslength < bitscutoff: # skip domains that are too short
			bitsRemovals += 1
			continue
		if evalue >= evaluecutoff: # skip domains with bad evalue
			evalueRemovals += 1
			continue

		# then count queries
		qseqid = lsplits[0]
		if donamechop: # for transdecoder peptides, |m.123 is needed for interval identification
			qseqid = qseqid.rsplit(donamechop,1)[0]
		querynamedict[qseqid] = True
		if swissprot:
		# blast outputs swissprot proteins as: sp|P0DI82|TPC2B_HUMAN
			sseqid = sseqid.split("|")[2] # should change to TPC2B_HUMAN
		else:
			sseqid = sseqid.replace("|","")
		hitDictCounter[sseqid] += 1

		backframe = False
		if hitstart > hitend: # for cases where transcript has backwards hit
			hitstart, hitend = hitend, hitstart # invert positions for calculation
			backframe = True # also change the strand
			backframecounts += 1

		genomeintervals = [] # to have empty iterable
		# convert protein positions to transcript nucleotide, as needed
		# protein position 1 becomes nucleotide position 1, position 2 becomes nucleotide 4, 3 to 7
		hitstart = (hitstart - 1) * multiplier + 1
		hitend = hitend * multiplier # end is necessarily the end of a codon
		scaffold = genescaffold.get(qseqid, None)
		strand = genestrand.get(qseqid, None)
		if backframe: # reassign strand if match is backwards
			strand = "+" if strand=="-" else "-"
		# convert transcript nucleotide to genomic nucleotide, and split at exon bounds
		if strand=='+':
			genomeintervals = get_intervals(geneintervals[qseqid], hitstart, hitlength, doreverse=False)
		elif strand=='-': # implies '-'
			genomeintervals = get_intervals(geneintervals[qseqid], hitstart, hitlength, doreverse=True)
		else: # strand is None
			print >> sys.stderr, "WARNING: cannot retrieve strand for {}".format(qseqid)
			continue
		intervalcounts += len(genomeintervals)
		if not len(genomeintervals):
			print >> sys.stderr, "WARNING: no intervals for {} in {}".format(sseqid, qseqid)
			intervalproblems += 1
			continue
		for interval in genomeintervals:
			# thus ID appears as qseqid.sseqid.number, so avic1234.avGFP.1, and uses ID in most browsers
			print >> sys.stdout, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t.\tID={7}.{8}.{9};Target={8} {10} {11} {12}".format(scaffold, programname, outputtype, interval[0], interval[1], bitscore, strand, qseqid, sseqid, hitDictCounter[sseqid], lsplits[8], lsplits[9], "-" if backframe else "+")
	print >> sys.stderr, "# Counted {} lines".format(linecounter), time.asctime()
	print >> sys.stderr, "# Removed {} hits by shortness".format(shortRemovals), time.asctime()
	print >> sys.stderr, "# Removed {} hits by bitscore".format(bitsRemovals), time.asctime()
	print >> sys.stderr, "# Removed {} hits by evalue".format(evalueRemovals), time.asctime()
	print >> sys.stderr, "# Found {} hits for {} queries".format(sum(hitDictCounter.values()), len(querynamedict) ), time.asctime()
	if backframecounts:
		print >> sys.stderr, "# {} hits are antisense".format(backframecounts), time.asctime()
	if intervalcounts:
		print >> sys.stderr, "# Wrote {} domain intervals".format(intervalcounts), time.asctime()
	if intervalproblems:
		print >> sys.stderr, "# {} genes have hits extending beyond gene bounds".format(intervalproblems), time.asctime()
	# NO RETURN

def get_intervals(intervals, domstart, domlength, doreverse=True):
	'''return a list of intervals with genomic positions for the feature'''
	# example domain arrangement for forward strand
	# intervals from     50,101 127,185 212,300
	# protein domain     71,101 127,185 212,256
	#      in nucleotides  31      59      45
	#      in amino acids  10.3    19.6    15 = 45
	# for domstart at 22 and domlength of 135
	# basestostart is always from transcript N-terminus nucleotide
	# so for forward transcripts, basestostart would be 22, so that 50+22-1=71
	basestostart = int(domstart) # this value always should be 1 or greater
	genomeintervals = [] # will contain a list of tuples
	for interval in sorted(intervals, key=lambda x: x[0], reverse=doreverse):
		intervallength = interval[1]-interval[0]+1 # corrected number of bases
		if basestostart >= intervallength: # ignore intervals before the start of the domain
		#	print >> sys.stderr, interval, domstart, basestostart, domlength, intervallength
			basestostart -= intervallength
		# in example, 101-50+1 = 52, 22 < 52, so else
		else: # bases to start is fewer than length of the interval, meaning domain must start here
			if doreverse: # reverse strand domains
				# if domain continues past an interval, domstart should be equal to interval[1]
				domstart = interval[1] - basestostart + 1 # correct for base numbering at end of interval
				if domstart - interval[0] + 1 >= domlength: # if the remaining part of the domain ends before the start of the interval
					# then define the last boundary and return the interval list
					genomebounds = (domstart-domlength+1, domstart) # subtract remaining length
					genomeintervals.append(genomebounds)
					return genomeintervals
				else:
					genomebounds = (interval[0], domstart)
					genomeintervals.append(genomebounds)
					domlength -= (domstart - interval[0] + 1)
					basestostart = 1 # start at the next interval
			else: # for forward stranded domains
				domstart = interval[0] + basestostart - 1 # correct for base numbering
				if interval[1] - domstart + 1 >= domlength: # if the remaining part of the domain ends before the end of the interval
					# then define the last boundary and return the interval list
					genomebounds = (domstart, domstart+domlength-1) # add remaining length for last interval
					genomeintervals.append(genomebounds)
					return genomeintervals
				else:
					genomebounds = (domstart, interval[1])
					genomeintervals.append(genomebounds)
					domlength -= (interval[1] - domstart + 1)
					basestostart = 1 # next domstart should be interval[0] for next interval
			if domlength < 1: # catch for if all domain length is accounted for
				return genomeintervals
	else:
		print >> sys.stderr, "WARNING: cannot finish protein at {} for {} in {}".format(domstart, domlength, intervals)
		return genomeintervals

def write_line(outlist, wayout):
	outline = "\t".join(outlist)
	print >> wayout, outline

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast', help="tabular blast results file")
	parser.add_argument('-d','--database', help="db proteins in fasta format")
	parser.add_argument('-g','--genes', help="query genes or proteins in gff format")
	parser.add_argument('-p','--program', help="blast program for 2nd column in output [BLASTX]", default="BLASTX")
	parser.add_argument('-t','--type', help="gff type or method [protein_match]", default="protein_match")
	parser.add_argument('-D','--delimiter', help="optional delimiter for protein names, cuts off end split")
	parser.add_argument('-c','--coverage-cutoff', type=float, help="query coverage cutoff for filtering [0.1]", default=0.1)
	parser.add_argument('-e','--evalue-cutoff', type=float, help="evalue cutoff [1e-3]", default=1e-3)
	parser.add_argument('-s','--score-cutoff', type=float, help="bitscore/length cutoff for filtering [0.1]", default=0.1)
	parser.add_argument('-F','--filter', action="store_true", help="filter low quality matches")
	parser.add_argument('-S','--swissprot', action="store_true", help="db sequences have swissprot headers")
	parser.add_argument('-T','--transdecoder', action="store_true", help="use presets for TransDecoder genome gff")
	parser.add_argument('-x','--exons', action="store_true", help="exons define coding sequence")
	parser.add_argument('-v','--verbose', action="store_true", help="extra output")
	args = parser.parse_args(argv)

	protlendb = make_seq_length_dict(args.database)
	geneintervals, genestrand, genescaffold =  gtf_to_intervals(args.genes, args.exons, args.transdecoder)

	parse_tabular_blast(args.blast, args.coverage_cutoff, args.evalue_cutoff, args.score_cutoff, args.program, args.type, args.delimiter, args.swissprot, protlendb, geneintervals, genestrand, genescaffold)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
