#!/usr/bin/env python
#
# pfam2gff.py v1.0 created 2016-03-03
#
# Sequence Ontology terms from:
# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo

'''
pfam2gff.py  last modified 2021-05-18

    EXAMPLE USAGE:
    to convert to protein gff, where domains are protein coordinates
pfam2gff.py -i proteins.pfam.tab > proteins.pfam.gtf

    to convert go genome gff, where domains are in genomic coordinates
    must use -g with transcript coordinates, can be exon or CDS
pfam2gff.py -i proteins.pfam.tab -g genes.gff > genome.pfam.gtf

    if using TransDecoder for peptide prediction, use mode -T
    and use the TransDecoder genome GFF file for -g
    this corrects for the CDS IDs as ID=cds.gene123 by removing the cds.
pfam2gff.py -i transdecoder.pfam.tab -g transdecoder_genome.gff3 -T > transdecoder.pfam.gff

    if using AUGUSTUS, use -T, and -d __ for extract_features.py proteins
pfam2gff.py -i augustus.prots.pfam.tab -g augustus.exons.gff -T -d __ > augustus.pfam.gff
    GENERATE AUGUSTUS EXONS BY:
reformatgff.py -a augustus.gff > augustus.exons.gff
    GENERATE AUGUSTUS PROTS BY:
extract_features.py -p augustus.prots.fasta augustus.gff

    GENERATE PFAM TABULAR BY:
hmmscan --cpu 4 --domtblout proteins.pfam.tab ~/PfamScan/data/Pfam-A.hmm transcripts.transdecoder.pep > proteins.pfam.log
'''

import os
import sys
import time
import argparse
import re
import gzip
from collections import defaultdict
from itertools import chain

def cds_to_intervals(gtffile, genesplit, keepexons, transdecoder, jgimode, nogenemode):
	'''convert protein or gene intervals from gff to dictionary where mrna IDs are keys and lists of intervals are values'''
	# this should probably be a class
	geneintervals = defaultdict(list)
	genestrand = {}
	genescaffold = {} 

	commentlines = 0
	linecounter = 0
	transcounter = 0
	exoncounter = 0
	# allow gzipped files
	if gtffile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing gff from {} as gzipped  ".format(gtffile) + time.asctime() + os.linesep)
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing gff from {}  ".format(gtffile) + time.asctime() + os.linesep)

	# extract gene or CDS information
	for line in opentype(gtffile,'r'):
		line = line.strip()
		if line: # ignore empty lines
			if line[0]=="#": # count comment lines, just in case
				commentlines += 1
			else:
				linecounter += 1
				lsplits = line.split("\t")
				scaffold = lsplits[0]
				feature = lsplits[2]
				strand = lsplits[6]
				attributes = lsplits[8]

				if attributes.find("ID")==0: # indicates gff3 format
					geneid = re.search('ID=([\w.|-]+)', attributes).group(1)
				elif attributes.find("Parent")==0: # gff3 format but no ID
					geneid = re.search('Parent=([\w.|-]+)', attributes).group(1)
				elif attributes.find("gene_id")==0: # indicates gtf format
					geneid = re.search('transcript_id "([\w.|-]+)";', attributes).group(1)
				elif jgimode and attributes.find("name")==0: # for JGI like: name "fgeneshTA2_pg.C_scaffold_1000001";
					geneid = re.search('name "([\w.|-]+)";', attributes).group(1)
				if transdecoder: # meaning CDS IDs will start with cds.gene.123|m.1
					geneid = geneid.replace("cds.","") # simply remove the cds.
					geneid = geneid.replace(".cds","") # also works for AUGUSTUS
				if genesplit: # if splitting, split all gene IDs
					geneid = geneid.split(genesplit,1)[0]
				if feature=="transcript" or feature=="mRNA":
					transcounter += 1
					genestrand[geneid] = strand
					genescaffold[geneid] = scaffold
				elif feature=="CDS" or (keepexons and feature=="exon"):
					exoncounter += 1
					boundaries = ( int(lsplits[3]), int(lsplits[4]) )
					if nogenemode: # gtf contains only exon and CDS, so get gene info from each CDS
						geneid = re.search('transcript_id "([\w.|-]+)";', attributes).group(1)
						# this may reassign multiple times
						genestrand[geneid] = strand
						genescaffold[geneid] = scaffold
					geneintervals[geneid].append(boundaries)
			#		sys.stderr.write("{} {} {}\n".format(geneid, scaffold, boundaries) )
	sys.stderr.write("# Counted {} lines and {} comments  ".format(linecounter, commentlines) + time.asctime() + os.linesep)
	sys.stderr.write("# Counted {} exons for {} transcripts  ".format(exoncounter, transcounter) + time.asctime() + os.linesep)
	sys.stderr.write("# Gene IDs taken as {} from {}\n".format(geneid, attributes) )
	return geneintervals, genestrand, genescaffold

def parse_pfam_domains(pfamtabular, evaluecutoff, lengthcutoff, programname, outputtype, donamechop, debugmode=False, jgimode=False, geneintervals=None, genestrand=None, genescaffold=None):
	'''parse domains from hmm domtblout and write to stdout as protein gff or genome gff'''
	domaincounter = 0
	protnamedict = {}
	evalueRemovals = 0
	shortRemovals = 0
	intervalproblems = 0
	intervalcounts = 0
	writecount = 0
	# for protein GFF, keep domains in dict for later sorting by position
	protboundstoline = defaultdict(dict)

	# allow gzipped files
	if pfamtabular.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing hmmscan PFAM tabular {} as gzipped  ".format(pfamtabular) + time.asctime() + os.linesep)
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing hmmscan PFAM tabular {}  ".format(pfamtabular) + time.asctime() + os.linesep)

	for line in opentype(pfamtabular, 'rt'):
		line = line.strip()
		if not line or line[0]=="#": # skip comment lines
			continue # also catch for empty line, which would cause IndexError
		domaincounter += 1
		lsplits = line.split(None, 22)
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
#      0                 1         2     3                    4         5        6       7     8     9  10   11       12        13     14    15     16   17     18   19     20  21  22
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
		targetname = lsplits[0]
		pfamacc = lsplits[1].rsplit('.',1)[0] # accs as PF00530.13, so chop off .13
		#queryid = lsplits[3].rsplit('|',1)[0] # chop transcript from transdecoder peptide
		queryid = lsplits[3] # for transdecoder peptides, |m.123 is needed for interval identification
		if donamechop:
			queryid = queryid.rsplit(donamechop,1)[0]
		protnamedict[queryid] = True
		evalue = float(lsplits[11]) # or [6] for full seq or [12]
		domscore = lsplits[13]
		domstart = int(lsplits[17])
		domend = int(lsplits[18])
		domnumber = lsplits[9]
		# include target description in the Name, change disallowed characters to -
		targetdescription = lsplits[22].replace("=","-").replace(",","-").replace(";","-")
		targetdescription = targetdescription.replace(" ","_")
		domainlength = domend - domstart + 1 # bases 1 to 6 should have length 6
		fractioncov = domainlength/float(lsplits[2])
		# filter matches by evalue, length, etc
		if fractioncov < lengthcutoff: # skip domains that are too short
			shortRemovals += 1
			continue
		if evalue >= evaluecutoff: # skip domains with bad evalue
			evalueRemovals += 1
			continue

		### FOR GENOME GFF ###
		if geneintervals: # if gene intervals are given in genomic coordinates
			genomeintervals = [] # to have empty iterable
			# convert domain protein positions to transcript nucleotide
			# protein position 1 becomes nucleotide position 1, position 2 becomes nucleotide 4, 3 to 7
			domstart = (domstart - 1) * 3 + 1
			domend = domend * 3 # end is necessarily the end of a codon
			domainlength_nucl = domend - domstart + 1 # recalculate for nucleotide coordinates
			scaffold = genescaffold.get(queryid, None)
			strand = genestrand.get(queryid, None)
			# convert transcript nucleotide to genomic nucleotide, and split at exon bounds
			if strand=='+':
				genomeintervals = get_intervals(geneintervals[queryid], domstart, domainlength_nucl, doreverse=False)
			elif strand=='-': # implies '-'
				genomeintervals = get_intervals(geneintervals[queryid], domstart, domainlength_nucl, doreverse=True)
			elif strand=='.': # strand is specified as '.'
				sys.stderr.write("WARNING: no strand given for {}, using forward\n".format(queryid) )
				genomeintervals = get_intervals(geneintervals[queryid], domstart, domainlength_nucl, doreverse=False)
			else: # strand is None, meaning queryid is not in genestrand dict
				sys.stderr.write("WARNING: no match in GFF for query {}\n".format(queryid) )
				continue
			intervalcounts += len(genomeintervals)
			if not len(genomeintervals):
				sys.stderr.write("WARNING: no intervals for {} in {}\n".format(targetname, queryid) )
				intervalproblems += 1
				continue

			writecount += 1

			# make Parent feature
			allpositions = list(chain(*genomeintervals))
			parentstart = min(allpositions)
			parentend = max(allpositions)
			# ID consists of: query gene, "target name" of up to 21 characters, domain number 
			# ID=g1.t1.VWA.1
			# Name consists of: PFAM accession, target name, target description
			# Name=PF00092.VWA.von_Willebrand_factor_type_A_domain
			outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t.\tID={10}.{8}.{9};Name={7}.{8}.{11}\n".format(scaffold, programname, outputtype, parentstart, parentend, domscore, strand, pfamacc, targetname, domnumber, queryid, targetdescription)
			sys.stdout.write(outline+os.linesep)
			# make child features for each interval
			for interval in genomeintervals:
				# thus ID appears as protein.targetname.number,
				# so avic1234.G2F.1, and uses ID in most browsers
				outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t.\tParent={10}.{8}.{9};Name={7}.{8}.{11}\n".format(scaffold, programname, outputtype, interval[0], interval[1], domscore, strand, pfamacc, targetname, domnumber, queryid, targetdescription)
				sys.stdout.write(outline)

		### FOR PROTEIN GFF ###
		else: # for protein GFF, make outline for later sorting
			writecount += 1
			boundaries = (domstart,domend)
			if debugmode:
				bitlength = float(lsplits[13])/domainlength
				outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{9:.3f}\t{10:.3f}\tID={6}.{7}.{8};Name={6}.{7}.{8}\n".format(queryid, programname, outputtype, domstart, domend, domscore, pfamacc, targetname, domnumber, bitlength, fractioncov)
			else:
				outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t.\t.\tID={6}.{7}.{8};Name={6}.{7}.{8}\n".format(queryid, programname, outputtype, domstart, domend, domscore, pfamacc, targetname, domnumber)
			protboundstoline[queryid][boundaries] = outline
	sys.stderr.write("# Found {} domains for {} proteins, wrote {}  {}\n".format(domaincounter, len(protnamedict), writecount, time.asctime() ) )
	if geneintervals: # in genome GFF mode, check if any CDS intervals were actually collected
		if intervalcounts:
			sys.stderr.write("# Wrote {} domain intervals\n".format(intervalcounts) )
		else:
			sys.stderr.write("WARNING: NO DOMAINS WRITTEN, CHECK OPTIONS -d AND -D\n")
	sys.stderr.write("# Removed {} domain hits by shortness\n".format(shortRemovals) )
	sys.stderr.write("# Removed {} domain hits by evalue\n".format(evalueRemovals) )
	if intervalproblems:
		sys.stderr.write("# {} genes have domains extending beyond gene bounds\n".format(intervalproblems) )
	if protboundstoline: # should be empty unless in protein GFF mode, meaning no genomic intervals
		for protid, boundlines in protboundstoline.items(): # sort proteins by start position
			for bounds in sorted(boundlines.keys()):
				sys.stdout.write(boundlines[bounds])
	# NO RETURN

def get_intervals(intervals, domstart, domlength, doreverse=True):
	'''return a list of intervals with genomic positions for the domain'''
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
		#	print(interval, domstart, basestostart, domlength, intervallength+os.linesep)
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
				# first interval, domstart should be 50+22-1 = 71
				# second interval, domstart should be 127+1-1 = 127
				# third interval, domstart should be 212+1-1 = 212
				domstart = interval[0] + basestostart - 1 # correct for base numbering
				if interval[1] - domstart + 1 >= domlength: # if the remaining part of the domain ends before the end of the interval
					# then define the last boundary and return the interval list
					genomebounds = (domstart, domstart+domlength-1) # add remaining length for last interval
					genomeintervals.append(genomebounds)
					return genomeintervals
				else:
					# first interval, genomebounds should be 71,101
					# domlength should be 135 - (101-71+1 = 31) = 104
					# second interval, genomebounds should be interval[0], domstart, so 127
					# domlength should be 104 - (185-127+1 = 59) = 45
					genomebounds = (domstart, interval[1])
					genomeintervals.append(genomebounds)
					domlength -= (interval[1] - domstart + 1)
					basestostart = 1 # next domstart should be interval[0] for next interval
			if domlength < 1: # catch for if all domain length is accounted for
				return genomeintervals
	else:
		sys.stderr.write("WARNING: cannot finish domain at {} for {} in {}\n".format(domstart, domlength, intervals) )
		return genomeintervals

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="PFAM domain information as hmmscan tabular")
	parser.add_argument('-d','--gene-delimiter', help="optional delimiter for gene names in GFF, cuts off end split")
	parser.add_argument('-D','--prot-delimiter', help="optional delimiter for protein names in PFAM table, cuts off end split")
	parser.add_argument('-e','--evalue', type=float, default=1e-1, help="evalue cutoff for domain filtering [1e-1]")
	parser.add_argument('-g','--genes', help="genes or proteins in gff format")
	parser.add_argument('-J','--JGI', action="store_true", help="use presets for JGI format files")
	parser.add_argument('-l','--length-cutoff', type=float, default=0.3, help="minimum length coverage of domains [0.3]")
	parser.add_argument('-n','--no-genes', action="store_true", help="genes are not defined, get gene ID for each exon")
	parser.add_argument('-p','--program', help="program for 2nd column in output [hmmscan]", default="hmmscan")
	parser.add_argument('-t','--type', help="gff type [PFAM]", default="PFAM")
	# this could more properly be SO:0000349 protein_match
	parser.add_argument('-T','--transdecoder', action="store_true", help="use presets for TransDecoder genome gff")
	parser.add_argument('-x','--exons', action="store_true", help="exons define coding sequence")
	parser.add_argument('--debug', action="store_true", help="debug some output options")
	args = parser.parse_args(argv)

	if args.genes:
		geneintervals, genestrand, genescaffold = cds_to_intervals(args.genes, args.gene_delimiter, args.exons, args.transdecoder, args.JGI, args.no_genes)
		parse_pfam_domains(args.input, args.evalue, args.length_cutoff, args.program, args.type, args.prot_delimiter, args.debug, args.JGI, geneintervals, genestrand, genescaffold)
	else: # assume protein gff
		parse_pfam_domains(args.input, args.evalue, args.length_cutoff, args.program, args.type, args.prot_delimiter, args.debug)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
